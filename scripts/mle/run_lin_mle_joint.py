#!/usr/bin/env python3
"""Fit a joint MAGeCK MLE model to the Lin OLA/DMSO3 screen.

Run beside lin_counts_OLA_DMSO3.count.txt:
    python run_lin_mle_joint.py
    python run_lin_mle_joint.py --dry-run

All 25 samples enter one model: three T0 samples, categorical effects for
T9/T12/T15/T18, and a shared OLA-versus-DMSO effect. The day indicators apply to
both treatment and control samples. T0 rows contain only the baseline term.

The drug coefficient is constant across days. The model has no treatment-by-day
interaction, replicate block or longitudinal covariance term. Sample counts
remain separate. Per-timepoint fits are archived in scripts/mle/previous/.

Whitespace and quote characters in guide/gene IDs receive reversible aliases
for the MAGeCK parser. Counts and guide membership are preserved; derived CSVs
restore the original gene IDs.
"""

import argparse
import csv
from datetime import datetime
from decimal import Decimal, InvalidOperation
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
TIMEPOINTS = {
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
NORM_METHOD = "median"
THREADS = 2
PERMUTATION_ROUNDS = 2  # MAGeCK default.
FDR_THRESHOLD = 0.05
EFFECT = "treatment_vs_control"
SCRIPT_VERSION = "2026-10-15-lin-joint-common-effect"


def make_design():
    """Categorical time effects plus a shared drug effect; no time trend."""
    days = list(TIMEPOINTS)
    header = ["Samples", "baseline", *[f"day_{day}" for day in days], EFFECT]
    rows = [[sample, 1, *([0] * len(days)), 0] for sample in BASELINE_SAMPLES]
    if not BASELINE_SAMPLES:
        raise ValueError("The configured Lin joint model requires the supplied T0 samples.")
    for index, day in enumerate(days):
        day_pattern = [int(i == index) for i in range(len(days))]
        for group, drug in [("control", 0), ("treatment", 1)]:
            if not TIMEPOINTS[day][group]:
                raise ValueError(f"Both DMSO and OLA samples are required for {day}.")
            rows.extend([sample, 1, *day_pattern, drug]
                        for sample in TIMEPOINTS[day][group])
    samples = [row[0] for row in rows]
    if len(samples) != len(set(samples)):
        raise ValueError("A sample appears more than once in the joint design.")
    return samples, header, rows


def integer_count(value):
    """Accept integer-valued exports such as 10.0; never round fractional counts."""
    try:
        number = Decimal(value.strip())
    except InvalidOperation as error:
        raise ValueError("A count is not numeric.") from error
    if not number.is_finite() or number < 0 or number != number.to_integral_value():
        raise ValueError("Nonnegative integer-valued raw counts are required.")
    return int(number)


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
                counts = [integer_count(value) for value in row[2:]]
            except ValueError as error:
                raise ValueError(f"Nonnegative raw integer counts required at line {line_no}; "
                                 "do not use normalized or log-transformed counts.") from error
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


def prepare_mageck_identifiers(rows):
    """Use reversible aliases instead of merging, trimming or dropping IDs."""
    mappings, records = [], []
    for column, kind in [(0, "sgRNA"), (1, "Gene")]:
        originals = {row[column] for row in rows}
        used = set(originals)
        mapping = {}
        serial = 0
        for name in sorted(originals):
            if not any(char.isspace() or char == '"' for char in name):
                continue
            while True:
                serial += 1
                alias = f"MAGECK_SAFE_{kind}_{serial:06d}"
                if alias not in used:
                    break
            used.add(alias)
            mapping[name] = alias
            records.append({"identifier_type": kind, "original_id": name, "mageck_id": alias})
        mappings.append(mapping)
    compatible = [[mappings[0].get(row[0], row[0]), mappings[1].get(row[1], row[1]),
                   *row[2:]] for row in rows]
    reverse_genes = {alias: original for original, alias in mappings[1].items()}
    return compatible, records, reverse_genes


def check_mageck_text_input(path, samples, expected_guides):
    """Check the tokenization used by MAGeCK 0.5.9.5 before launching a fit."""
    seen = 0
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.strip().split()
            if line_no == 1:
                if fields[2:] != samples or len(fields) != len(samples) + 2:
                    raise ValueError(f"MAGeCK cannot parse the sample header: {path}")
                continue
            if len(fields) != len(samples) + 2:
                raise ValueError(f"MAGeCK-incompatible fields at {path}, line {line_no}.")
            try:
                for value in fields[2:]:
                    int(value)
            except ValueError as error:
                raise ValueError(f"MAGeCK-incompatible count at {path}, line {line_no}.") from error
            seen += 1
    if seen != expected_guides:
        raise ValueError(f"Guide-count mismatch after preparing {path}: {seen} != {expected_guides}")


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


def summarize(prefix, reverse_genes):
    source = Path(str(prefix) + ".gene_summary.txt")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, results = reader.fieldnames or [], list(reader)
    required = ["Gene", f"{EFFECT}|beta", f"{EFFECT}|p-value", f"{EFFECT}|fdr"]
    if set(required) - set(fields) or not results:
        raise ValueError("MAGeCK did not produce the expected joint treatment-effect results.")
    if reverse_genes:
        for row in results:
            row["Gene"] = reverse_genes.get(row["Gene"], row["Gene"])
        write_table(Path(str(prefix) + ".gene_summary.original_ids.txt"), fields,
                    ([row[field] for field in fields] for row in results))
    for row in results:
        beta = result_number(row[f"{EFFECT}|beta"])
        fdr = result_number(row[f"{EFFECT}|fdr"])
        row["delta_beta"] = beta  # Direct coefficient from this joint fit.
        row["positive_selection"] = math.isfinite(beta) and beta > 0
        row["significant_positive"] = (row["positive_selection"] and math.isfinite(fdr)
                                       and 0 <= fdr < FDR_THRESHOLD)
    results.sort(key=lambda row: row["delta_beta"] if math.isfinite(row["delta_beta"])
                 else -math.inf, reverse=True)
    for rank, row in enumerate(results, 1):
        row["positive_rank"] = rank
    extra = ["delta_beta", "positive_rank", "positive_selection", "significant_positive"]
    hits = [row for row in results if row["significant_positive"]]
    tag = str(FDR_THRESHOLD).replace(".", "")
    write_records(Path(str(prefix) + "_ranked.csv"), fields + extra, results)
    write_records(Path(str(prefix) + f"_positive_FDR{tag}.csv"), fields + extra, hits)
    print(f"\nJoint model: {len(results)} genes/entries; {len(hits)} with common "
          f"OLA beta > 0 and permutation FDR < {FDR_THRESHOLD}", flush=True)
    print("Gene\tcommon_OLA_beta\tpermutation_FDR", flush=True)
    for row in results[:10]:
        print(f"{row['Gene']}\t{row['delta_beta']:.5g}\t{row[f'{EFFECT}|fdr']}", flush=True)
    return {"genes": len(results), "positive_hits": len(hits),
            "effect_column": f"{EFFECT}|beta", "fdr_column": f"{EFFECT}|fdr",
            "fdr_threshold": FDR_THRESHOLD}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count-table", type=Path, default=COUNT_FILE)
    parser.add_argument("--dry-run", action="store_true",
                        help="Prepare the single joint design and counts without fitting.")
    args = parser.parse_args()
    count_file = args.count_table.expanduser().resolve()
    if not count_file.is_file():
        raise FileNotFoundError(f"Count table not found: {count_file}")
    samples, design_header, design_rows = make_design()
    header, rows, qc = read_counts(count_file, set(samples))
    compatible_rows, identifier_map, reverse_genes = prepare_mageck_identifiers(rows)
    mageck = shutil.which("mageck")
    if not args.dry_run and not mageck:
        raise RuntimeError("MAGeCK not found. Run: conda activate mageck_env")
    version = "not checked (dry run)"
    if not args.dry_run:
        check = subprocess.run([mageck, "--version"], check=True, text=True, capture_output=True)
        version = (check.stdout + check.stderr).strip()

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    suffix = "_dry_run" if args.dry_run else ""
    out = FOLDER / "MAGeCK_MLE_results" / f"LIN_JOINT_{stamp}{suffix}"
    out.mkdir(parents=True, exist_ok=False)
    selected = out / "selected.count.txt"
    indices = [0, 1, *[header.index(sample) for sample in samples]]
    write_table(selected, [header[i] for i in indices],
                ([row[i] for i in indices] for row in compatible_rows))
    check_mageck_text_input(selected, samples, len(rows))
    design_file = out / "design.txt"
    write_table(design_file, design_header, design_rows)
    write_records(out / "count_qc.csv", list(qc[0]), qc)
    write_records(out / "identifier_map.csv", ["identifier_type", "original_id", "mageck_id"],
                  identifier_map)
    prefix = out / "lin_joint_OLA_vs_DMSO3"
    command = [mageck or "mageck", "mle", "-k", str(selected), "-d", str(design_file),
               "--norm-method", NORM_METHOD, "--threads", str(THREADS),
               "--permutation-round", str(PERMUTATION_ROUNDS), "-n", str(prefix)]
    metadata = {
        "script_version": SCRIPT_VERSION, "input": str(count_file),
        "input_sha256": hashlib.sha256(count_file.read_bytes()).hexdigest(),
        "model": "one common OLA-versus-DMSO effect adjusted for categorical day",
        "baseline_samples": BASELINE_SAMPLES, "timepoint_samples": TIMEPOINTS,
        "selected_samples_in_order": samples, "design_columns": design_header[1:],
        "assumptions": [
            "One common treatment effect across T9, T12, T15 and T18; no treatment-by-day interactions.",
            "Categorical time effects are unrestricted across the four days.",
            "No replicate-block or longitudinal residual-correlation model is included.",
            "A/B/C pairing is not inferred from column labels; shared lineages need experimental metadata.",
        ],
        "normalization": NORM_METHOD, "threads": THREADS,
        "permutation_rounds": PERMUTATION_ROUNDS, "mageck_version": version,
        "python_version": sys.version, "dry_run": args.dry_run,
        "target_effect": f"{EFFECT}|beta", "target_fdr": f"{EFFECT}|fdr",
        "fdr_threshold": FDR_THRESHOLD, "identifier_map": str(out / "identifier_map.csv"),
        "aliased_identifiers": len(identifier_map), "command": command,
    }
    (out / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print("Script version:", SCRIPT_VERSION, flush=True)
    print(f"Validated {len(rows)} guides; ONE joint fit with {len(samples)} sample columns.", flush=True)
    print("T0: 3; DMSO: 11; OLA: 11. Target: common OLA effect adjusted for day.", flush=True)
    print("Output directory:", out, flush=True)
    print(f"Identifier aliases: {len(identifier_map)}; all {len(rows)} guides retained.", flush=True)
    print("\nDesign:\n" + "\t".join(design_header), flush=True)
    for row in design_rows:
        print("\t".join(map(str, row)), flush=True)
    print("\nCommand:", shlex.join(command), flush=True)
    if args.dry_run:
        print("Dry run complete. One joint input prepared; MAGeCK was not run.", flush=True)
        return
    subprocess.run(command, check=True)
    summary = summarize(prefix, reverse_genes)
    (out / "joint_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("\nJoint MLE completed. Results:", out, flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

