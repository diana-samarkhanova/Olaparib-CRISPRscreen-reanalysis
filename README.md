# Olaparib CRISPR screen reanalysis

Code and dataset records for *Beyond Homologous Recombination: Functional Modules of Olaparib Resistance Revealed by CRISPR Screens*.

The repository contains RRA commands for six of the ten publications in the review and MLE workflows for Lau, Lin, Tuano and Clements. Input count tables, complete run logs and study results are not included. [Data availability and verification status](docs/data_status.md) describes the current coverage.

## Repository contents

| Path | Contents |
| --- | --- |
| `scripts/mageck/` | MAGeCK count and RRA commands by study |
| `scripts/mle/` | MLE workflows and the archived Lin per-timepoint configuration |
| `scripts/collect_hits.py` | RRA hit selection and recurrence counting |
| `metadata/studies.tsv` | Studies and input sources |
| `metadata/comparisons.tsv` | RRA comparisons and recurrence-unit decisions |
| `metadata/mle_comparisons.tsv` | MLE samples and model parameters |
| `metadata/mle_designs/` | MLE design matrices |
| `provenance/` | Original command records and script checksums |
| `docs/` | Methods, dataset records and data availability |

## RRA

FASTQ inputs are processed with `mageck count`, followed by `mageck test`. Existing sgRNA count tables enter directly at `mageck test`.

Activate the MAGeCK environment and arrange the input files as shown in the scripts. Each shell script takes the study working directory as its only argument:

```bash
# FASTQ input
bash scripts/mageck/lau_count.sh /path/to/1LAU
bash scripts/mageck/lau_test.sh /path/to/1LAU

# Existing count table
bash scripts/mageck/tuano_test.sh /path/to/3TUANO
```

The scripts retain the sample selections, arguments and output prefixes recorded in `provenance/command_records.json`. Lin uses the final recorded count command. Existing outputs are protected from overwriting. Dataset-specific settings and unresolved records are listed in [dataset records](docs/analysis_notes.md).

## MLE

[MLE analysis](docs/mle_analysis.md) describes installation, input files, model definitions and outputs. Lin uses one shared OLA effect adjusted for sampling day. Tuano cell models and Clements genotypes are fitted separately.

Positive MLE hits have `treatment_vs_control|beta > 0` and permutation `treatment_vs_control|fdr < 0.05`. MLE recurrence aggregation is not implemented in this repository.

## RRA recurrence analysis

`collect_hits.py` reads MAGeCK RRA gene-summary files and counts supporting recurrence units. Multiple comparisons assigned to the same `recurrence_unit_id` contribute once per gene.

Each included comparison in `metadata/comparisons.tsv` requires its result path, model ID, `include_in_recurrence=yes`, `review_status=verified`, recurrence-unit ID and decision note. Excluded comparisons require `include_in_recurrence=no` and a decision note. The script rejects unresolved entries.

```bash
python3 scripts/collect_hits.py \
  --manifest metadata/comparisons.tsv \
  --data-root /path/to/study-directories \
  --outdir results/run_01 \
  --fdr 0.05 --min-screens 2 \
  --lfc-policy report-only
```

This example selects by `pos|fdr` and reports `pos|lfc`. The alternative `--lfc-policy require-positive` also requires `pos|lfc > 0`. The LFC rule used for the manuscript gene list is not documented in the archived records.

Outputs include per-comparison statistics, hit counts, recurrent genes, supporting studies and the manifest used. Python 3.9 or later is required.

## Reactome

[Reactome analysis](docs/reactome_workflow.md) describes the enrichment metadata and module-assignment tables. Seven modules were based on enrichment results; the eighth group was informed by the literature. The original enrichment export and analysis settings are not archived.

## Tests

```bash
python3 -m unittest discover -s tests -v
```

The tests use synthetic data to check hit selection, recurrence counting and input validation.
