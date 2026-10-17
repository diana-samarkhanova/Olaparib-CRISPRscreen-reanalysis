#!/usr/bin/env python3
"""Run five separate MAGeCK MLE comparisons for the Tuano OL/DMSO screens.

Place this script beside tuano_count_edited.txt, activate mageck_env,
and run: python run_tuano_mle.py
Optional: python run_tuano_mle.py --dry-run
Optional: python run_tuano_mle.py --comparisons HMLE B80_T17
Optional: python run_tuano_mle.py --count-table /path/to/count_table.txt

No separate T0 samples are present in the supplied header. For each comparison:
    Samples   baseline  treatment_vs_control
    DMSO_1    1         0
    DMSO_2    1         0
    OL_1      1         1
    OL_2      1         1

DMSO is the statistical reference, not an assumed day-zero measurement.
Each cell model is fitted separately. Sample columns remain
separate, without averaging or a replicate-block term. Labels 1/2 alone are not
taken as evidence of a matched-pair design. Positive beta means relative
enrichment under OL versus DMSO within that comparison. There is no separately
estimated DMSO-versus-T0 effect because this table has no T0 samples.
FDR is within each comparison, not jointly over all five tests. B80-T5 and
B80-T17 are distinct cell lines; define the recurrence unit before integration.
Special characters in sample names are preserved; output folder names are safe.

Only Python's standard library and an existing MAGeCK installation are needed.
This script validates and selects raw counts; it does not pre-normalize counts.
MAGeCK 0.5.9.5 splits TSV input on all whitespace. IDs containing whitespace
or double quotes therefore receive collision-free temporary aliases. All
counts, guides and gene groupings are preserved; identifier_map.csv records
the mapping. Ranked CSV files use the original gene names. When genes were
aliased, an additional .gene_summary.original_ids.txt is written; the original
MAGeCK .gene_summary.txt retains the exact engine output with temporary IDs.
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
COUNT_FILE_NAMES = ("tuano_count_edited.txt", "tuano_count_edited.tsv", "tuano_count_edited")
COMPARISONS = {
    "HMLE": {
        "cell_model": "HMLE", "timepoint": None,
        "control": ["HMLE_DMSO_1", "HMLE_DMSO_2"],
        "treatment": ["HMLE_OL_1", "HMLE_OL_2"],
    },
    "mesHMLE": {
        "cell_model": "mesHMLE", "timepoint": None,
        "control": ["mesHMLE_DMSO_1", "mesHMLE_DMSO_2"],
        "treatment": ["mesHMLE_OL_1", "mesHMLE_OL_2"],
    },
    "B80_T5": {
        "cell_model": "B80-T5", "timepoint": None,
        "control": ["B80_T5_DMSO_1", "B80_T5_DMSO_2"],
        "treatment": ["B80_T5_OL_1", "B80_T5_OL_2"],
    },
    "B80_T17": {
        "cell_model": "B80-T17", "timepoint": None,
        "control": ["B80_T17_DMSO_1", "B80_T17_DMSO_2"],
        "treatment": ["B80_T17_OL_1", "B80_T17_OL_2"],
    },
    "K5_K19": {
        "cell_model": "K5+/K19-", "timepoint": None,
        "control": ["K5+/K19-_DMSO_1", "K5+/K19-_DMSO_2"],
        "treatment": ["K5+/K19-_OL_1", "K5+/K19-_OL_2"],
    },
}
NORM_METHOD = "median"
THREADS = 2
PERMUTATION_ROUNDS = 2  # MAGeCK default, as in the previous template-based run.
FDR_THRESHOLD = 0.05
EFFECT = "treatment_vs_control"
SCRIPT_VERSION = "2026-10-17-median-checked"


def require_median_counts(qc):
    """Reject the automatic total fallback used by MAGeCK 0.5.9.5."""
    if NORM_METHOD != "median":
        raise ValueError("This workflow requires NORM_METHOD = 'median'.")
    blocked = [f"{item['sample']} ({item['zero_percent']:.2f}% zeros)"
               for item in qc if item.get("selected", True)
               and 100 * item["zero_count_guides"] > 45 * item["guides"]]
    if blocked:
        raise ValueError(
            "Median normalization cannot be accepted: MAGeCK 0.5.9.5 switches "
            "to total when any selected sample has >45% zero-count guides. "
            + "; ".join(blocked) + ". No MLE fit was run. Review the count table "
            "and normalization choice; do not remove zero rows just to force median."
        )


def verify_median_log(prefix):
    """Check the engine log before exporting hits; retain evidence per fit."""
    log_path = Path(str(prefix) + ".log")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    fallback = "Switch to total read count normalization." in text
    factors = [line for line in text.splitlines() if "Final size factor:" in line]
    verified = bool(factors) and not fallback
    report = {"requested": NORM_METHOD, "log": str(log_path),
              "median_verified_from_log": verified, "total_fallback": fallback,
              "final_size_factor_log_lines": factors}
    Path(str(prefix) + ".normalization.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    if not verified:
        raise RuntimeError(f"Median normalization was not verified in {log_path}. "
                           "This run must not be used as a median-normalized result.")


def make_design(comparison):
    groups = [COMPARISONS[comparison]["control"], COMPARISONS[comparison]["treatment"]]
    if any(not group for group in groups):
        raise ValueError("Both control and treatment groups are required.")
    samples = [sample for group in groups for sample in group]
    if len(samples) != len(set(samples)):
        raise ValueError(f"A sample belongs to multiple groups in {comparison}.")
    header = ["Samples", "baseline", EFFECT]
    patterns = [[1, 0], [1, 1]]
    rows = [[sample, *pattern] for group, pattern in zip(groups, patterns)
            for sample in group]
    return samples, header, rows


def find_count_file(explicit_path=None):
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Count table does not exist: {path}")
        return path
    candidates = [FOLDER / name for name in COUNT_FILE_NAMES if (FOLDER / name).is_file()]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise ValueError("Multiple count tables found. Select one using --count-table PATH: "
                         + ", ".join(str(path) for path in candidates))
    raise FileNotFoundError(f"Place this script beside tuano_count_edited.txt in {FOLDER}, "
                            "or specify --count-table PATH.")


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


def summarize(prefix, comparison, reverse_genes=None):
    with Path(str(prefix) + ".gene_summary.txt").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, results = reader.fieldnames or [], list(reader)
    required = ["Gene", f"{EFFECT}|beta", f"{EFFECT}|p-value", f"{EFFECT}|fdr"]
    if set(required) - set(fields):
        raise ValueError(f"Unexpected MLE output columns. Required: {required}")
    if not results:
        raise ValueError(f"No gene results were produced for {comparison}.")
    if reverse_genes:
        for row in results:
            row["Gene"] = reverse_genes.get(row["Gene"], row["Gene"])
        write_table(Path(str(prefix) + ".gene_summary.original_ids.txt"), fields,
                    ([row[field] for field in fields] for row in results))
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
    print(f"\n{comparison}: {len(results)} genes/entries; "
          f"{len(hits)} with beta > 0 and permutation FDR < {FDR_THRESHOLD}", flush=True)
    print("Gene\tbeta\tpermutation_FDR", flush=True)
    for row in results[:10]:
        print(f"{row['Gene']}\t{row['delta_beta']:.5g}\t{row[f'{EFFECT}|fdr']}", flush=True)
    group = COMPARISONS[comparison]
    return {"comparison": comparison, "cell_model": group["cell_model"],
            "timepoint": group["timepoint"], "genes": len(results), "positive_hits": len(hits),
            "fdr_threshold": FDR_THRESHOLD, "fdr_column": f"{EFFECT}|fdr",
            "result_directory": str(prefix.parent)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate and prepare all inputs without fitting MAGeCK.")
    parser.add_argument("--comparisons", nargs="+", choices=list(COMPARISONS),
                        default=list(COMPARISONS), help="Default: all five comparisons")
    parser.add_argument("--count-table", type=Path,
                        help="Default: tuano_count_edited.txt/.tsv beside this script.")
    args = parser.parse_args()
    comparisons = list(dict.fromkeys(args.comparisons))
    designs = {name: make_design(name) for name in comparisons}
    required_samples = {sample for samples, _, _ in designs.values() for sample in samples}
    count_file = find_count_file(args.count_table)
    header, rows, qc = read_counts(count_file, required_samples)
    compatible_rows, identifier_map, reverse_genes = prepare_mageck_identifiers(rows)
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
    out = FOLDER / "MAGeCK_MLE_results" / f"TUANO_{stamp}{suffix}"
    out.mkdir(parents=True, exist_ok=False)
    write_records(out / "count_qc.csv", list(qc[0]), qc)
    write_records(out / "identifier_map.csv", ["identifier_type", "original_id", "mageck_id"],
                  identifier_map)
    metadata = {
        "script_version": SCRIPT_VERSION,
        "input": str(count_file), "input_sha256": hashlib.sha256(count_file.read_bytes()).hexdigest(),
        "separate_T0_samples": [], "reference_group": "DMSO within each comparison",
        "comparisons": {name: COMPARISONS[name] for name in comparisons},
        "normalization": NORM_METHOD, "threads": THREADS,
        "permutation_rounds": PERMUTATION_ROUNDS, "fdr_column": f"{EFFECT}|fdr",
        "fdr_threshold": FDR_THRESHOLD, "mageck_version": version,
        "python_version": sys.version, "dry_run": args.dry_run,
        "identifier_map": str(out / "identifier_map.csv"),
        "aliased_identifiers": len(identifier_map),
        "identifier_handling": "Whitespace/double-quote IDs receive unique reversible aliases. "
                               "No rows are dropped or merged. Counts are unchanged. "
                               "Ranked CSVs and gene_summary.original_ids.txt use original gene IDs; "
                               "raw MAGeCK outputs retain engine IDs.",
        "design_note": "Separate two-group models: DMSO is reference; no T0 samples. "
                       "No replicate-block effect. No sample averaging. "
                       "B80-T5 and B80-T17 are distinct cell lines; recurrence unit requires a separate decision. "
                       "FDR is within each comparison; no across-comparison adjustment.",
        "commands": {},
    }
    jobs = []
    for comparison, (samples, design_header, design_rows) in designs.items():
        folder = out / comparison
        folder.mkdir()
        indices = [0, 1, *[header.index(sample) for sample in samples]]
        selected = folder / "selected.count.txt"
        write_table(selected, [header[i] for i in indices],
                    ([row[i] for i in indices] for row in compatible_rows))
        check_mageck_text_input(selected, samples, len(rows))
        design_file = folder / "design.txt"
        write_table(design_file, design_header, design_rows)
        prefix = folder / f"tuano_{comparison}_OL_vs_DMSO"
        command = [mageck or "mageck", "mle", "-k", str(selected),
                   "-d", str(design_file), "--norm-method", NORM_METHOD,
                   "--threads", str(THREADS), "--permutation-round", str(PERMUTATION_ROUNDS),
                   "-n", str(prefix)]
        (folder / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
        metadata["commands"][comparison] = command
        jobs.append((comparison, prefix, command))
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Validated {len(rows)} guides across {len(header) - 2} sample columns.", flush=True)
    print("Script version:", SCRIPT_VERSION, flush=True)
    print("Output directory:", out, flush=True)
    print(f"MAGeCK-compatible identifier aliases: {len(identifier_map)}; "
          f"all {len(rows)} guides retained.", flush=True)
    for entry in identifier_map[:10]:
        print(f"{entry['identifier_type']}: {entry['original_id']!r} -> {entry['mageck_id']}", flush=True)
    for comparison in comparisons:
        group = COMPARISONS[comparison]
        print(f"{comparison}: DMSO={len(group['control'])}, OL={len(group['treatment'])}", flush=True)
    for item in qc:
        if item["selected"]:
            print(f"{item['sample']}: {item['total_counts']} counts; "
                  f"{item['zero_percent']:.2f}% zero-count guides", flush=True)
    require_median_counts(qc)
    if args.dry_run:
        print("Dry run complete. All inputs prepared; MAGeCK was not run.", flush=True)
        return

    summaries = []
    for comparison, prefix, command in jobs:
        print("\nRunning", comparison, ":", shlex.join(command), flush=True)
        subprocess.run(command, check=True)
        verify_median_log(prefix)
        summaries.append(summarize(prefix, comparison, reverse_genes))
    write_records(out / "comparison_summary.csv", list(summaries[0]), summaries)
    print("\nAll requested comparisons completed. Results:", out, flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
