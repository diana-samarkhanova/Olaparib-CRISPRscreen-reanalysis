#!/usr/bin/env python3
"""Fit MAGeCK MLE for Lau OL_1 versus Control.

Run beside lau_counts.count.txt:
    python run_lau_mle.py
    python run_lau_mle.py --dry-run

The configured comparison has one sample per group and no separate T0.
OL_2 exclusion and the identity of Control are unresolved in the source records.
Optional BASELINE_SAMPLES adds a separate baseline group; group membership must
be unique. Outputs include the design, selected raw counts, QC and ranked hits.
"""

import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


# Configuration
FOLDER = Path(__file__).resolve().parent
COUNT_FILES = ["lau_counts.count.txt"]
BASELINE_SAMPLES = []  # No separate T0 column.
CONTROL_SAMPLES = ["Control"]
TREATMENT_SAMPLES = ["OL_1"]
NORM_METHOD = "total"  # Matches the recorded Lau RRA normalization fallback.
THREADS = 2
PERMUTATION_ROUNDS = 2  # MAGeCK default.
FDR_THRESHOLD = 0.05
EFFECT = "treatment_vs_control"


def make_design():
    groups = [BASELINE_SAMPLES, CONTROL_SAMPLES, TREATMENT_SAMPLES]
    if not CONTROL_SAMPLES or not TREATMENT_SAMPLES:
        raise ValueError("Both control and treatment samples are required.")
    samples = [sample for group in groups for sample in group]
    if len(samples) != len(set(samples)):
        raise ValueError("Sample names must be unique across groups.")
    if BASELINE_SAMPLES:
        header = ["Samples", "baseline", "control", EFFECT]
        patterns = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
    else:
        header = ["Samples", "baseline", EFFECT]
        patterns = [[], [1, 0], [1, 1]]
    rows = [[sample, *pattern] for group, pattern in zip(groups, patterns)
            for sample in group]
    return samples, header, rows


def prepare_counts(path, selected_samples, output_path):
    """Check raw counts and preserve every guide while selecting columns."""
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader, [])
        if len(header) < 4 or len(header) != len(set(header)):
            raise ValueError("Expected unique headers: sgRNA, Gene, samples.")
        if header[0].lower() != "sgrna" or header[1].lower() != "gene":
            raise ValueError("The first two columns must be sgRNA and Gene.")
        missing = sorted(set(selected_samples) - set(header[2:]))
        if missing:
            raise ValueError(f"Missing sample columns: {missing}")
        indices = [header.index(sample) for sample in selected_samples]
        totals, zeros = [0] * (len(header) - 2), [0] * (len(header) - 2)
        seen, n_guides = set(), 0
        with output_path.open("x", encoding="utf-8", newline="") as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow([*header[:2], *selected_samples])
            for line_no, row in enumerate(reader, 2):
                if not row or not any(row):
                    continue
                if len(row) != len(header) or not row[0] or not row[1]:
                    raise ValueError(f"Invalid count-table row at line {line_no}.")
                if row[0] in seen:
                    raise ValueError(f"Duplicate sgRNA ID at line {line_no}: {row[0]}")
                seen.add(row[0])
                values = [float(value) for value in row[2:]]
                if any(not math.isfinite(value) or value < 0 or not value.is_integer()
                       for value in values):
                    raise ValueError(f"Expected nonnegative integer counts at line {line_no}.")
                for i, value in enumerate(values):
                    totals[i] += int(value)
                    zeros[i] += value == 0
                writer.writerow([*row[:2], *[int(values[i - 2]) for i in indices]])
                n_guides += 1
    if not n_guides:
        raise ValueError("The count table contains no guides.")
    qc = [{"sample": sample, "selected": sample in selected_samples,
           "guides": n_guides, "total_counts": total,
           "zero_count_guides": zero, "zero_percent": round(100 * zero / n_guides, 4)}
          for sample, total, zero in zip(header[2:], totals, zeros)]
    if any(item["selected"] and item["total_counts"] == 0 for item in qc):
        raise ValueError("A selected sample has no counts.")
    return qc


def write_rows(path, fields, rows):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize(prefix):
    with Path(str(prefix) + ".gene_summary.txt").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, results = reader.fieldnames or [], list(reader)
    required = ["Gene", f"{EFFECT}|beta", f"{EFFECT}|p-value", f"{EFFECT}|fdr"]
    if set(required) - set(fields):
        raise ValueError(f"Unexpected MLE result columns. Required: {required}")
    for row in results:
        beta, fdr = float(row[f"{EFFECT}|beta"]), float(row[f"{EFFECT}|fdr"])
        # Directly fitted treatment-versus-control coefficient.
        row["delta_beta"] = beta
        row["positive_selection"] = math.isfinite(beta) and beta > 0
        row["significant_positive"] = (row["positive_selection"]
                                       and math.isfinite(fdr)
                                       and 0 <= fdr < FDR_THRESHOLD)
    results.sort(key=lambda row: row["delta_beta"]
                 if math.isfinite(row["delta_beta"]) else -math.inf, reverse=True)
    for rank, row in enumerate(results, 1):
        row["positive_rank"] = rank
    extra = ["delta_beta", "positive_rank", "positive_selection", "significant_positive"]
    hits = [row for row in results if row["significant_positive"]]
    write_rows(Path(str(prefix) + "_ranked.csv"), fields + extra, results)
    write_rows(Path(str(prefix) + "_positive_FDR005.csv"), fields + extra, hits)
    print(f"\nGenes/entries: {len(results)}; beta > 0 and permutation FDR < 0.05: {len(hits)}")
    print("Gene\tbeta\tpermutation_FDR")
    for row in results[:10]:
        print(f"{row['Gene']}\t{row['delta_beta']:.5g}\t{row[f'{EFFECT}|fdr']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    samples, design_header, design_rows = make_design()
    mageck = shutil.which("mageck")
    if not args.dry_run and not mageck:
        raise RuntimeError("MAGeCK is not on PATH. Run: conda activate mageck_env")
    version = "not checked (dry run)"
    if not args.dry_run:
        check = subprocess.run([mageck, "--version"], check=True, text=True,
                               capture_output=True)
        version = (check.stdout + check.stderr).strip()
    if len(CONTROL_SAMPLES) < 2 or len(TREATMENT_SAMPLES) < 2:
        print("NOTE: fewer than two sample columns in at least one comparison group.")
        print("MLE FDR values are model-based and do not establish biological replication.")

    for filename in COUNT_FILES:
        count_file = FOLDER / filename
        if not count_file.is_file():
            raise FileNotFoundError(f"Put this script beside the count table: {count_file}")
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        out = FOLDER / "MAGeCK_MLE_results" / f"{count_file.stem}_{stamp}"
        out.mkdir(parents=True, exist_ok=False)
        selected = out / "selected.count.txt"
        qc = prepare_counts(count_file, samples, selected)
        write_rows(out / "count_qc.csv", list(qc[0]), qc)
        for item in qc:
            status = "selected" if item["selected"] else "not selected"
            print(f"{item['sample']}: {status}; {item['total_counts']} counts; "
                  f"{item['zero_percent']:.2f}% zero guides")
        design_file = out / "design.txt"
        with design_file.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(design_header)
            writer.writerows(design_rows)
        prefix = out / "mle_treatment_vs_control"
        command = [mageck or "mageck", "mle", "-k", str(selected),
                   "-d", str(design_file), "--norm-method", NORM_METHOD,
                   "--threads", str(THREADS), "--permutation-round", str(PERMUTATION_ROUNDS),
                   "-n", str(prefix)]
        metadata = {"input": str(count_file),
                    "input_sha256": hashlib.sha256(count_file.read_bytes()).hexdigest(),
                    "baseline_samples": BASELINE_SAMPLES, "control_samples": CONTROL_SAMPLES,
                    "treatment_samples": TREATMENT_SAMPLES, "normalization": NORM_METHOD,
                    "mageck_version": version, "python_version": sys.version,
                    "threads": THREADS, "permutation_rounds": PERMUTATION_ROUNDS,
                    "fdr_column": f"{EFFECT}|fdr", "fdr_threshold": FDR_THRESHOLD,
                    "dry_run": args.dry_run, "command": command,
                    "sample_selection_note": "OL_1 vs Control follows the recorded Lau RRA. "
                    "OL_2 exclusion rationale and count QC require review before final inference."}
        (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                                encoding="utf-8")
        (out / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
        print("\nDesign:")
        print("\t".join(design_header))
        for row in design_rows:
            print("\t".join(map(str, row)))
        print("\nCommand:", shlex.join(command))
        if args.dry_run:
            print("\nDry run finished. Inputs saved to:", out)
            continue
        subprocess.run(command, check=True)
        summarize(prefix)
        print("\nResults saved to:", out)


if __name__ == "__main__":
    main()
