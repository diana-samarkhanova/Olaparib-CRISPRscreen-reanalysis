#!/usr/bin/env python3
"""Fit separate MAGeCK MLE models for Clements HeLa WT and BRCA2KO.

Run beside Clements_HeLa_WT.count.txt and Clements_HeLa_BRCA2KO.count.txt:
    python run_clements_mle.py
    python run_clements_mle.py --dry-run
    python run_clements_mle.py --datasets HeLa_BRCA2KO
    python run_clements_mle.py --data-dir /path/to/counts
    python run_clements_mle.py --permutation-rounds 10

Each genotype has two DMSO and two OL samples, with DMSO as the reference.
The fits have no separate T0, replicate-block or genotype-by-drug interaction
term. Positive beta indicates enrichment under OL. FDR is adjusted within each
fit. Defaults: median normalization, two threads and two permutation rounds.

Source: Clements et al. 2020, DOI 10.1038/s41467-020-19961-w.
MOESM4 supplies HeLa WT and MOESM5 supplies BRCA2KO. On the worksheet
"sgRNA scores and read counts", A/B map to sgRNA/Gene and N/P/O/Q map to
DMSO_1/DMSO_2/OL_1/OL_2. Each export retains all 77,441 guides and raw counts.
Reference hashes identify the recorded exports; input validation checks file
structure and count values, not the experimental identity of the samples.

Both inputs are validated before fitting. Integer-valued Excel exports are
accepted; fractional, missing, negative and nonfinite counts are rejected.
Reversible aliases handle whitespace/quotes in guide and gene IDs. Derived CSVs
restore original IDs; maps and engine outputs are retained for each dataset.
Every run writes a new timestamped results directory.

Documentation: https://sourceforge.net/p/mageck/wiki/usage/
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
DATASETS = {
    "HeLa_WT": {
        "file": "Clements_HeLa_WT.count.txt",
        "genotype": "WT",
        "source_supplement": "41467_2020_19961_MOESM4_ESM.xlsx",
        "reference_export_sha256":
            "df6efb09c7326f88a57b7772ec8a5dcbf551715da4f125106100d1c2b55df5d1",
    },
    "HeLa_BRCA2KO": {
        "file": "Clements_HeLa_BRCA2KO.count.txt",
        "genotype": "BRCA2KO",
        "source_supplement": "41467_2020_19961_MOESM5_ESM.xlsx",
        "reference_export_sha256":
            "35a1197d760a301143ce46953c31b2f081db6a19a6f018629a51d41b6b980b12",
    },
}
CONTROL = ["DMSO_1", "DMSO_2"]
TREATMENT = ["OL_1", "OL_2"]
SAMPLES = CONTROL + TREATMENT
EFFECT = "treatment_vs_control"
NORM_METHOD = "median"
THREADS = 2
PERMUTATION_ROUNDS = 2
FDR_THRESHOLD = 0.05
SCRIPT_VERSION = "2026-10-15-clements-v1"


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("Use an integer greater than zero.")
    return number


def make_design():
    header = ["Samples", "baseline", EFFECT]
    rows = [[sample, 1, 0] for sample in CONTROL]
    rows += [[sample, 1, 1] for sample in TREATMENT]
    return header, rows


def integer_count(value):
    """Parse integer-valued exports exactly; never round, replace or impute."""
    try:
        number = Decimal(value.strip())
    except InvalidOperation as error:
        raise ValueError("not a numeric count") from error
    if not number.is_finite():
        raise ValueError("missing or nonfinite count")
    if number < 0:
        raise ValueError("negative count")
    if number != number.to_integral_value():
        raise ValueError("fractional count")
    return int(number)


def read_counts(path):
    """Read strict TSV; sample identity, not column position, defines groups."""
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            header = next(reader, [])
            if len(header) != 6 or len(header) != len(set(header)):
                raise ValueError(
                    f"{path.name}: expected 6 unique TAB-separated columns: "
                    "sgRNA, Gene, DMSO_1, DMSO_2, OL_1, OL_2. "
                    "Do not supply gene-level scores or extra Excel columns."
                )
            if [value.lower() for value in header[:2]] != ["sgrna", "gene"]:
                raise ValueError(f"{path.name}: the first columns must be sgRNA and Gene.")
            if set(header[2:]) != set(SAMPLES):
                raise ValueError(
                    f"{path.name}: sample labels must be {SAMPLES}; found {header[2:]}. "
                    "Check spaces in the header and the original sample mapping."
                )
            totals, zeros = [0] * 4, [0] * 4
            seen, rows = set(), []
            for row in reader:
                line_no = reader.line_num
                if not row or not any(row):
                    continue
                if len(row) != len(header):
                    raise ValueError(
                        f"{path.name}, line {line_no}: {len(row)} fields instead of 6."
                    )
                for column in (0, 1):
                    if not row[column].strip():
                        raise ValueError(
                            f"{path.name}, line {line_no}: empty {header[column]}."
                        )
                if row[0] in seen:
                    raise ValueError(
                        f"{path.name}, line {line_no}: duplicate sgRNA {row[0]!r}."
                    )
                seen.add(row[0])
                counts = []
                for sample, token in zip(header[2:], row[2:]):
                    try:
                        counts.append(integer_count(token))
                    except ValueError as error:
                        raise ValueError(
                            f"{path.name}, line {line_no}, sample {sample}: "
                            f"value {token!r} is invalid ({error}). "
                            "Use original nonnegative integer read counts; "
                            "do not round normalized values or fill missing counts."
                        ) from error
                for i, value in enumerate(counts):
                    totals[i] += value
                    zeros[i] += value == 0
                rows.append([*row[:2], *counts])
    except UnicodeError as error:
        raise ValueError(
            f"{path.name}: cannot read UTF-8 text. Re-export as UTF-8 tab-separated "
            "text; do not rename an Excel workbook to .txt."
        ) from error
    except csv.Error as error:
        raise ValueError(f"{path.name}: malformed tab-separated text ({error}).") from error
    if not rows:
        raise ValueError(f"{path.name}: the table contains no sgRNAs.")
    qc = [
        {"sample": sample, "guides": len(rows), "total_counts": total,
         "zero_count_guides": zero, "zero_percent": round(100 * zero / len(rows), 4)}
        for sample, total, zero in zip(header[2:], totals, zeros)
    ]
    empty = [entry["sample"] for entry in qc if entry["total_counts"] == 0]
    if empty:
        raise ValueError(f"{path.name}: samples with no counts: {empty}.")
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

def summarize(prefix, dataset, reverse_genes):
    """Keep the complete fit plus both directions of significant selection."""
    path = Path(str(prefix) + ".gene_summary.txt")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, results = reader.fieldnames or [], list(reader)
    required = ["Gene", f"{EFFECT}|beta", f"{EFFECT}|p-value", f"{EFFECT}|fdr"]
    if set(required) - set(fields):
        raise ValueError(f"{path}: missing MLE result columns; required: {required}")
    if not results:
        raise ValueError(f"{dataset}: MAGeCK produced no gene results.")
    if reverse_genes:
        for row in results:
            row["Gene"] = reverse_genes.get(row["Gene"], row["Gene"])
        write_table(Path(str(prefix) + ".gene_summary.original_ids.txt"), fields,
                    ([row[field] for field in fields] for row in results))

    for row in results:
        beta = result_number(row[f"{EFFECT}|beta"])
        fdr = result_number(row[f"{EFFECT}|fdr"])
        # Directly fitted OL-versus-DMSO coefficient.
        row["delta_beta"] = beta
        row["positive_selection"] = math.isfinite(beta) and beta > 0
        row["negative_selection"] = math.isfinite(beta) and beta < 0
        significant = math.isfinite(fdr) and 0 <= fdr < FDR_THRESHOLD
        row["significant_positive"] = row["positive_selection"] and significant
        row["significant_negative"] = row["negative_selection"] and significant
    results.sort(key=lambda row: row["delta_beta"]
                 if math.isfinite(row["delta_beta"]) else -math.inf, reverse=True)
    for rank, row in enumerate(results, 1):
        row["positive_rank"] = rank
    positive = [row for row in results if row["significant_positive"]]
    negative = sorted((row for row in results if row["significant_negative"]),
                      key=lambda row: row["delta_beta"])
    extra = ["delta_beta", "positive_rank", "positive_selection", "negative_selection",
             "significant_positive", "significant_negative"]
    tag = str(FDR_THRESHOLD).replace(".", "")
    write_records(Path(str(prefix) + "_ranked.csv"), fields + extra, results)
    write_records(Path(str(prefix) + f"_positive_FDR{tag}.csv"), fields + extra, positive)
    write_records(Path(str(prefix) + f"_negative_FDR{tag}.csv"), fields + extra, negative)
    nonfinite = sum(not math.isfinite(row["delta_beta"]) for row in results)
    print(f"\n{dataset}: {len(results)} gene entries; "
          f"{len(positive)} positive and {len(negative)} negative hits "
          f"at permutation FDR < {FDR_THRESHOLD}.", flush=True)
    if nonfinite:
        print(f"WARNING: {nonfinite} nonfinite beta values retained in the full result.",
              flush=True)
    return {
        "dataset": dataset, "cell_line": "HeLa", "genotype": DATASETS[dataset]["genotype"],
        "contrast": "OL_vs_DMSO", "gene_entries": len(results),
        "positive_hits": len(positive), "negative_hits": len(negative),
        "nonfinite_betas": nonfinite, "fdr_column": f"{EFFECT}|fdr",
        "fdr_threshold": FDR_THRESHOLD, "result_directory": str(prefix.parent),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="Check both tables and prepare inputs without running MAGeCK.")
    parser.add_argument("--datasets", nargs="+", choices=list(DATASETS),
                        default=list(DATASETS), help="Default: both genotypes.")
    parser.add_argument("--data-dir", type=Path, default=FOLDER,
                        help="Directory containing both count tables; default: beside script.")
    parser.add_argument("--threads", type=positive_int, default=THREADS)
    parser.add_argument("--permutation-rounds", type=positive_int,
                        default=PERMUTATION_ROUNDS)
    args = parser.parse_args()
    names = list(dict.fromkeys(args.datasets))
    data_dir = args.data_dir.expanduser().resolve()
    print("Script version:", SCRIPT_VERSION, flush=True)
    prepared = {}

    # Validate all inputs before fitting.
    for dataset in names:
        settings = DATASETS[dataset]
        path = data_dir / settings["file"]
        if not path.is_file():
            raise FileNotFoundError(
                f"Count table not found: {path}\n"
                "Place this script in /mnt/e/DIANA/6CLEMENTS/counts, "
                "or specify --data-dir /path/to/counts."
            )
        header, rows, qc = read_counts(path)
        compatible, mapping, reverse_genes = prepare_mageck_identifiers(rows)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        matches = digest == settings["reference_export_sha256"]
        prepared[dataset] = {
            "path": path, "header": header, "rows": compatible, "qc": qc,
            "mapping": mapping, "reverse_genes": reverse_genes,
            "sha256": digest, "matches_reference": matches,
            "gene_entries": len({row[1] for row in rows}),
            "all_zero_guides": sum(all(value == 0 for value in row[2:]) for row in rows),
        }
        print(f"{dataset}: validated {len(rows)} unique guides, "
              f"{prepared[dataset]['gene_entries']} gene labels, 4 samples.", flush=True)
        if matches:
            print("  Exact match to the previously validated published-count export.",
                  flush=True)
        else:
            print("  NOTE: file differs from the saved reference export. Format/count "
                  "checks passed; verify any changes against the source table.", flush=True)
        print(f"  Temporary identifier aliases: {len(mapping)}; "
              f"all-zero guides retained: {prepared[dataset]['all_zero_guides']}.",
              flush=True)
        for item in qc:
            print(f"  {item['sample']}: {item['total_counts']} counts; "
                  f"{item['zero_percent']:.2f}% zero-count guides.", flush=True)

    mageck = shutil.which("mageck")
    if not args.dry_run and not mageck:
        raise RuntimeError("MAGeCK not found. Run: conda activate mageck_env")
    version = "not checked (dry run)"
    if not args.dry_run:
        check = subprocess.run([mageck, "--version"], check=True, capture_output=True,
                               text=True)
        version = (check.stdout + check.stderr).strip()

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    suffix = "_dry_run" if args.dry_run else ""
    out = data_dir / "MAGeCK_MLE_results" / f"CLEMENTS_{stamp}{suffix}"
    out.mkdir(parents=True, exist_ok=False)
    print("Output directory:", out, flush=True)
    print(f"Normalization: {NORM_METHOD}; threads: {args.threads}; "
          f"permutation rounds: {args.permutation_rounds}.", flush=True)
    metadata = {
        "script_version": SCRIPT_VERSION, "created_local": datetime.now().isoformat(),
        "source_doi": "10.1038/s41467-020-19961-w",
        "source_worksheet": "sgRNA scores and read counts",
        "source_column_mapping": {
            "sgRNA": "A", "Gene": "B", "DMSO_1": "N", "DMSO_2": "P",
            "OL_1": "O", "OL_2": "Q",
        },
        "conversion_note": "Published individual-sample read counts N:Q; "
                           "statistic columns C:M not used. No prior normalization, "
                           "averaging or filtering in the recorded conversion.",
        "reference_export_guides_per_file": 77441,
        "normalization": NORM_METHOD, "threads": args.threads,
        "permutation_rounds": args.permutation_rounds,
        "fdr_column": f"{EFFECT}|fdr", "fdr_threshold": FDR_THRESHOLD,
        "mageck_version": version, "python_version": sys.version, "dry_run": args.dry_run,
        "design_note": "Separate OL-vs-DMSO model per genotype. DMSO is reference. "
                       "No T0 samples, sample averaging or inferred replicate blocks. "
                       "No genotype-by-treatment interaction test. "
                       "FDR is within each fit, without across-dataset adjustment.",
        "identifier_note": "Reversible temporary aliases for whitespace/quote IDs. "
                           "No sgRNAs are removed or merged. Counts are unchanged. "
                           "Derived CSVs restore gene IDs; engine output stays intact.",
        "datasets": {}, "commands": {},
    }
    jobs = []
    all_qc = []
    design_header, design_rows = make_design()
    for dataset, info in prepared.items():
        folder = out / dataset
        folder.mkdir()
        indices = [0, 1, *[info["header"].index(sample) for sample in SAMPLES]]
        selected = folder / "selected.count.txt"
        write_table(selected, ["sgRNA", "Gene", *SAMPLES],
                    ([row[i] for i in indices] for row in info["rows"]))
        check_mageck_text_input(selected, SAMPLES, len(info["rows"]))
        design_file = folder / "design.txt"
        write_table(design_file, design_header, design_rows)
        write_records(folder / "identifier_map.csv",
                      ["identifier_type", "original_id", "mageck_id"], info["mapping"])
        all_qc.extend({"dataset": dataset, **item} for item in info["qc"])
        prefix = folder / f"clements_{dataset}_OL_vs_DMSO"
        command = [
            mageck or "mageck", "mle", "-k", str(selected), "-d", str(design_file),
            "--norm-method", NORM_METHOD, "--threads", str(args.threads),
            "--permutation-round", str(args.permutation_rounds), "-n", str(prefix),
        ]
        (folder / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
        metadata["commands"][dataset] = command
        metadata["datasets"][dataset] = {
            **DATASETS[dataset], "input_path": str(info["path"]),
            "input_sha256": info["sha256"], "matches_reference_export": info["matches_reference"],
            "guides": len(info["rows"]), "gene_entries": info["gene_entries"],
            "all_zero_guides_retained": info["all_zero_guides"],
            "identifier_aliases": len(info["mapping"]),
            "samples": SAMPLES,
        }
        jobs.append((dataset, prefix, command, info["reverse_genes"]))
    write_records(out / "count_qc.csv", list(all_qc[0]), all_qc)
    (out / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    # Both selected inputs passed the MAGeCK parser check.
    del prepared
    if args.dry_run:
        print("Dry run complete. Inputs passed validation; MAGeCK was not run.",
              flush=True)
        return

    summaries = []
    summary_path = out / "comparison_summary.csv"
    for dataset, prefix, command, reverse_genes in jobs:
        print("\nRunning", dataset, ":", shlex.join(command), flush=True)
        subprocess.run(command, check=True)
        summaries.append(summarize(prefix, dataset, reverse_genes))
        # Save completed work even if a later MAGeCK fit fails.
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
            writer.writeheader()
            writer.writerows(summaries)
    print("\nAll requested comparisons completed. Results:", out, flush=True)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted. Completed files remain in the printed output directory.",
              file=sys.stderr)
        sys.exit(130)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

