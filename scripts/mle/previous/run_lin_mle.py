#!/usr/bin/env python3
"""Run four separate MAGeCK MLE comparisons for the Lin OLA/DMSO3 screen.

Place this script beside lin_counts_OLA_DMSO3.count.txt, activate mageck_env,
and run: python run_lin_mle.py
Optional: python run_lin_mle.py --dry-run
Optional: python run_lin_mle.py --timepoints T18

At each time point the model follows Code_mle.docx:
    Samples       baseline  control  treatment_vs_control
    T0            1         0        0
    same-day DMSO 1         1        0
    same-day OLA  1         1        1

All three T0 samples form a shared baseline in each separate analysis.
Sample columns are kept separate, without averaging or a replicate-block term.
The A/B/C labels alone are not taken as evidence of a matched-pair design.
Different time points are not treated as independent biological replicates.
Positive beta means relative enrichment under OLA versus same-day DMSO.
FDR is calculated within each comparison, not jointly over all four tests.

Only Python's standard library and an existing MAGeCK installation are needed.
This script validates and selects raw counts; it does not pre-normalize counts.
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


FOLDER = Path(__file__).resolve().parent
COUNT_FILE = FOLDER / "lin_counts_OLA_DMSO3.count.txt"
BASELINE_SAMPLES = ["T0_Screen3_C", "T0_Screen3_B", "T0_Screen3_A"]
COMPARISONS = {
    "T9": {
        "control": ["T9_DMSO3_B", "T9_DMSO3_A"],
        "treatment": ["T9_OLA_B", "T9_OLA_A"],
    },
    "T12": {
        "control": ["T12_DMSO3_C", "T12_DMSO3_B", "T12_DMSO3_A"],
        "treatment": ["T12_OLA_C", "T12_OLA_B", "T12_OLA_A"],
    },
    "T15": {
        "control": ["T15_DMSO3_C", "T15_DMSO3_B", "T15_DMSO3_A"],
        "treatment": ["T15_OLA_C", "T15_OLA_B", "T15_OLA_A"],
    },
    "T18": {
        "control": ["T18_DMSO3_C", "T18_DMSO3_B", "T18_DMSO3_A"],
        "treatment": ["T18_OLA_C", "T18_OLA_B", "T18_OLA_A"],
    },
}
NORM_METHOD = "median"  # Template/default; assess control-guide normalization if needed.
THREADS = 2
PERMUTATION_ROUNDS = 2  # MAGeCK default, as in the previous template-based run.
FDR_THRESHOLD = 0.05
EFFECT = "treatment_vs_control"


def make_design(timepoint):
    groups = [BASELINE_SAMPLES, COMPARISONS[timepoint]["control"],
              COMPARISONS[timepoint]["treatment"]]
    if any(not group for group in groups):
        raise ValueError("Separate T0, control and treatment groups are required.")
    samples = [sample for group in groups for sample in group]
    if len(samples) != len(set(samples)):
        raise ValueError(f"A sample belongs to multiple groups in {timepoint}.")
    header = ["Samples", "baseline", "control", EFFECT]
    patterns = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
    rows = [[sample, *pattern] for group, pattern in zip(groups, patterns)
            for sample in group]
    return samples, header, rows


def read_counts(path, required_samples):
    """Validate raw integer counts and preserve every guide, including zeros."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) < 4 or len(header) != len(set(header)):
            raise ValueError("Expected unique headers: sgRNA, Gene, samples.")
        if header[0].lower() != "sgrna" or header[1].lower() != "gene":
            raise ValueError("The first two count-table columns must be sgRNA and Gene.")
        missing = sorted(set(required_samples) - set(header[2:]))
        if missing:
            raise ValueError(f"Missing sample columns: {missing}")
        totals, zeros = [0] * (len(header) - 2), [0] * (len(header) - 2)
        seen, rows = set(), []
        for line_no, row in enumerate(reader, 2):
            if not row or not any(row):
                continue
            if len(row) != len(header) or not row[0].strip() or not row[1].strip():
                raise ValueError(f"Invalid count-table row at line {line_no}.")
            if row[0] in seen:
                raise ValueError(f"Duplicate sgRNA at line {line_no}: {row[0]}")
            seen.add(row[0])
            try:
                counts = [int(value) for value in row[2:]]
            except ValueError as error:
                raise ValueError(f"Raw integer counts required at line {line_no}.") from error
            if any(value < 0 for value in counts):
                raise ValueError(f"Negative count at line {line_no}.")
            for i, value in enumerate(counts):
                totals[i] += value
                zeros[i] += value == 0
            rows.append([*row[:2], *counts])
    if not rows:
        raise ValueError("The count table contains no guides.")
    qc = [{"sample": sample, "selected": sample in required_samples,
           "guides": len(rows), "total_counts": total, "zero_count_guides": zero,
           "zero_percent": round(100 * zero / len(rows), 4)}
          for sample, total, zero in zip(header[2:], totals, zeros)]
    if any(item["selected"] and item["total_counts"] == 0 for item in qc):
        raise ValueError("A selected sample has no counts.")
    return header, rows, qc


def write_table(path, header, rows, delimiter="\t"):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_records(path, fields, rows):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def result_number(value):
    if value is None or str(value).strip().lower() in {"", "na", "n/a", "none"}:
        return math.nan
    return float(value)


def summarize(prefix, timepoint):
    with Path(str(prefix) + ".gene_summary.txt").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, results = reader.fieldnames or [], list(reader)
    required = ["Gene", f"{EFFECT}|beta", f"{EFFECT}|p-value", f"{EFFECT}|fdr"]
    if set(required) - set(fields):
        raise ValueError(f"Unexpected MLE output columns. Required: {required}")
    if not results:
        raise ValueError(f"No gene results were produced for {timepoint}.")
    for row in results:
        beta = result_number(row[f"{EFFECT}|beta"])
        fdr = result_number(row[f"{EFFECT}|fdr"])
        # Direct fitted contrast, not a subtraction of separately fitted betas.
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
    tag = str(FDR_THRESHOLD).replace(".", "")
    write_records(Path(str(prefix) + "_ranked.csv"), fields + extra, results)
    write_records(Path(str(prefix) + f"_positive_FDR{tag}.csv"), fields + extra, hits)
    print(f"\n{timepoint}: {len(results)} genes/entries; "
          f"{len(hits)} with beta > 0 and permutation FDR < {FDR_THRESHOLD}", flush=True)
    print("Gene\tbeta\tpermutation_FDR", flush=True)
    for row in results[:10]:
        print(f"{row['Gene']}\t{row['delta_beta']:.5g}\t{row[f'{EFFECT}|fdr']}", flush=True)
    return {"timepoint": timepoint, "genes": len(results), "positive_hits": len(hits),
            "fdr_threshold": FDR_THRESHOLD, "fdr_column": f"{EFFECT}|fdr",
            "result_directory": str(prefix.parent)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate and prepare all inputs without fitting MAGeCK.")
    parser.add_argument("--timepoints", nargs="+", choices=list(COMPARISONS),
                        default=list(COMPARISONS), help="Default: T9 T12 T15 T18")
    args = parser.parse_args()
    timepoints = list(dict.fromkeys(args.timepoints))
    designs = {tp: make_design(tp) for tp in timepoints}
    required_samples = {sample for samples, _, _ in designs.values() for sample in samples}
    if not COUNT_FILE.is_file():
        raise FileNotFoundError(f"Place this script beside the count table: {COUNT_FILE}")
    header, rows, qc = read_counts(COUNT_FILE, required_samples)
    mageck = shutil.which("mageck")
    if not args.dry_run and not mageck:
        raise RuntimeError("MAGeCK not found. Run: conda activate mageck_env")
    version = "not checked (dry run)"
    if not args.dry_run:
        check = subprocess.run([mageck, "--version"], check=True, text=True,
                               capture_output=True)
        version = (check.stdout + check.stderr).strip()

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    suffix = "_dry_run" if args.dry_run else ""
    out = FOLDER / "MAGeCK_MLE_results" / f"LIN_{stamp}{suffix}"
    out.mkdir(parents=True, exist_ok=False)
    write_records(out / "count_qc.csv", list(qc[0]), qc)
    metadata = {
        "input": str(COUNT_FILE), "input_sha256": hashlib.sha256(COUNT_FILE.read_bytes()).hexdigest(),
        "baseline_samples": BASELINE_SAMPLES,
        "comparisons": {tp: COMPARISONS[tp] for tp in timepoints},
        "normalization": NORM_METHOD, "threads": THREADS,
        "permutation_rounds": PERMUTATION_ROUNDS, "fdr_column": f"{EFFECT}|fdr",
        "fdr_threshold": FDR_THRESHOLD, "mageck_version": version,
        "python_version": sys.version, "dry_run": args.dry_run,
        "design_note": "Separate per-timepoint models with all three T0 samples. "
                       "No replicate-block effect. No sample averaging. "
                       "The four timepoints are not independent screens. "
                       "FDR is within each comparison; no across-timepoint adjustment.",
        "commands": {},
    }
    jobs = []
    for timepoint, (samples, design_header, design_rows) in designs.items():
        folder = out / timepoint
        folder.mkdir()
        indices = [0, 1, *[header.index(sample) for sample in samples]]
        selected = folder / "selected.count.txt"
        write_table(selected, [header[i] for i in indices],
                    ([row[i] for i in indices] for row in rows))
        design_file = folder / "design.txt"
        write_table(design_file, design_header, design_rows)
        prefix = folder / f"lin_{timepoint}_OLA_vs_DMSO3"
        command = [mageck or "mageck", "mle", "-k", str(selected),
                   "-d", str(design_file), "--norm-method", NORM_METHOD,
                   "--threads", str(THREADS), "--permutation-round", str(PERMUTATION_ROUNDS),
                   "-n", str(prefix)]
        (folder / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
        metadata["commands"][timepoint] = command
        jobs.append((timepoint, prefix, command))
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Validated {len(rows)} guides across {len(header) - 2} sample columns.", flush=True)
    print("Output directory:", out, flush=True)
    for timepoint in timepoints:
        group = COMPARISONS[timepoint]
        print(f"{timepoint}: T0={len(BASELINE_SAMPLES)}, DMSO={len(group['control'])}, "
              f"OLA={len(group['treatment'])}", flush=True)
    if args.dry_run:
        print("Dry run complete. All inputs prepared; MAGeCK was not run.", flush=True)
        return

    summaries = []
    for timepoint, prefix, command in jobs:
        print("\nRunning", timepoint, ":", shlex.join(command), flush=True)
        subprocess.run(command, check=True)
        summaries.append(summarize(prefix, timepoint))
    write_records(out / "comparison_summary.csv", list(summaries[0]), summaries)
    print("\nAll requested comparisons completed. Results:", out, flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
