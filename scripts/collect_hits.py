#!/usr/bin/env python3
"""Select positive MAGeCK hits and count support across defined screen units."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate header fields: {path}")
        rows = list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError(f"Malformed TSV row: {path}")
        return reader.fieldnames, rows


def write_table(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value, name, context):
    try:
        result = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid {name} in {context}: {value!r}") from None
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {name} in {context}: {value!r}")
    return result


def collect(manifest, data_root, outdir, fdr, min_screens, lfc_policy):
    if not 0 < fdr <= 1:
        raise ValueError("FDR threshold must be greater than 0 and at most 1")
    if min_screens < 1:
        raise ValueError("Minimum recurrence must be at least 1")
    if lfc_policy not in {"report-only", "require-positive"}:
        raise ValueError("An explicit, supported LFC policy is required")

    manifest = Path(manifest)
    data_root = Path(data_root)
    outdir = Path(outdir)
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise ValueError("Output directory must be new or empty; preserve previous runs")
    fields, rows = read_table(manifest)
    required = {"comparison_id", "study_id", "model_id", "gene_summary_file",
                "include_in_recurrence", "recurrence_unit_id", "review_status", "decision_note"}
    if not required.issubset(fields):
        raise ValueError(f"Missing manifest fields: {sorted(required - set(fields))}")
    if not rows:
        raise ValueError("The comparison manifest has no rows")

    included, excluded = [], []
    comparison_ids = set()
    for row in rows:
        cid = row["comparison_id"].strip()
        if not cid or cid in comparison_ids:
            raise ValueError(f"Missing or duplicate comparison ID: {cid!r}")
        comparison_ids.add(cid)
        decision = row["include_in_recurrence"].strip().lower()
        if decision not in {"yes", "no"}:
            raise ValueError(f"Pending inclusion decision for {cid}")
        if not row["decision_note"].strip():
            raise ValueError(f"A decision note is required for {cid}")
        if decision == "no":
            excluded.append(row)
            continue
        if row["review_status"].strip() != "verified":
            raise ValueError(f"Included comparison is not verified: {cid}")
        for name in ("study_id", "model_id", "gene_summary_file", "recurrence_unit_id"):
            if not row[name].strip():
                raise ValueError(f"Missing {name} for included comparison {cid}")
        included.append(row)
    if not included:
        raise ValueError("No comparisons were explicitly included")

    evidence = []
    comparison_summary = []
    support = defaultdict(lambda: {"units": set(), "studies": set(), "comparisons": set(), "models": set()})
    used_paths = set()
    for comparison in included:
        cid = comparison["comparison_id"]
        summary_path = Path(comparison["gene_summary_file"])
        if not summary_path.is_absolute():
            summary_path = data_root / summary_path
        resolved = summary_path.resolve()
        if resolved in used_paths:
            raise ValueError(f"The same gene-summary file is assigned to multiple included comparisons: {summary_path}")
        used_paths.add(resolved)
        columns, genes = read_table(summary_path)
        required_columns = {"id", "pos|fdr", "pos|lfc"}
        if not required_columns.issubset(columns):
            raise ValueError(f"Missing positive-selection fields in {summary_path}: {sorted(required_columns - set(columns))}")
        if not genes:
            raise ValueError(f"Empty gene summary: {summary_path}")
        seen_genes = set()
        n_positive_fdr = n_hits = n_nonpositive_lfc = 0
        for gene_row in genes:
            gene = gene_row["id"].strip()
            if not gene or gene in seen_genes:
                raise ValueError(f"Missing or duplicate gene ID in {cid}: {gene!r}")
            seen_genes.add(gene)
            context = f"{cid}/{gene}"
            positive_fdr = numeric(gene_row["pos|fdr"], "pos|fdr", context)
            positive_lfc = numeric(gene_row["pos|lfc"], "pos|lfc", context)
            if not 0 <= positive_fdr <= 1:
                raise ValueError(f"FDR outside [0,1] in {context}")
            significant = positive_fdr < fdr
            hit = significant and (lfc_policy == "report-only" or positive_lfc > 0)
            n_positive_fdr += int(significant)
            n_nonpositive_lfc += int(significant and positive_lfc <= 0)
            n_hits += int(hit)
            evidence.append({
                "gene_id": gene, "study_id": comparison["study_id"],
                "model_id": comparison["model_id"], "comparison_id": cid,
                "recurrence_unit_id": comparison["recurrence_unit_id"],
                "pos_fdr": positive_fdr, "pos_lfc": positive_lfc,
                "positive_fdr_pass": int(significant), "hit": int(hit),
                "source_gene_summary": comparison["gene_summary_file"],
            })
            if hit:
                support[gene]["units"].add(comparison["recurrence_unit_id"])
                support[gene]["studies"].add(comparison["study_id"])
                support[gene]["comparisons"].add(cid)
                support[gene]["models"].add(comparison["model_id"])
        comparison_summary.append({
            "comparison_id": cid, "study_id": comparison["study_id"],
            "n_tested_genes": len(genes), "n_positive_fdr": n_positive_fdr,
            "n_hits_after_lfc_policy": n_hits,
            "n_nonpositive_lfc_among_positive_fdr": n_nonpositive_lfc,
        })

    recurrent = []
    for gene, hits in sorted(support.items()):
        if len(hits["units"]) < min_screens:
            continue
        recurrent.append({
            "gene_id": gene, "n_recurrence_units": len(hits["units"]),
            "n_publications": len(hits["studies"]), "n_supporting_comparisons": len(hits["comparisons"]),
            "recurrence_unit_ids": ";".join(sorted(hits["units"])),
            "study_ids": ";".join(sorted(hits["studies"])),
            "model_ids": ";".join(sorted(hits["models"])),
            "comparison_ids": ";".join(sorted(hits["comparisons"])),
        })
    # All data and scientific-decision checks run before any output is written.
    outdir.mkdir(parents=True, exist_ok=True)
    write_table(outdir / "comparison_gene_evidence.tsv", list(evidence[0]), evidence)
    write_table(outdir / "comparison_summary.tsv", list(comparison_summary[0]), comparison_summary)
    recurrent_fields = ["gene_id", "n_recurrence_units", "n_publications", "n_supporting_comparisons",
                        "recurrence_unit_ids", "study_ids", "model_ids", "comparison_ids"]
    write_table(outdir / "recurrent_genes.tsv", recurrent_fields, recurrent)
    (outdir / "recurrent_genes.txt").write_text("".join(row["gene_id"] + "\n" for row in recurrent), encoding="utf-8")
    write_table(outdir / "manifest_used.tsv", fields, rows)
    report = {
        "manifest": str(manifest), "fdr_threshold": fdr, "fdr_comparison": "strictly_less_than",
        "lfc_policy": lfc_policy, "minimum_recurrence_units": min_screens,
        "n_included_comparisons": len(included), "n_excluded_comparisons": len(excluded),
        "n_included_publications": len({x["study_id"] for x in included}),
        "n_included_recurrence_units": len({x["recurrence_unit_id"] for x in included}),
        "n_recurrent_genes": len(recurrent),
        "note": "Recurrence units were supplied in the verified manifest; independence is not inferred by this script.",
    }
    (outdir / "run_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fdr", type=float, required=True)
    parser.add_argument("--min-screens", type=int, required=True)
    parser.add_argument("--lfc-policy", choices=["report-only", "require-positive"], required=True)
    args = parser.parse_args()
    try:
        result = collect(args.manifest, args.data_root, args.outdir, args.fdr, args.min_screens, args.lfc_policy)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
