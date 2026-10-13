# Olaparib CRISPR screen reanalysis

This repository is being prepared to document the reanalysis underlying
**Beyond Homologous Recombination: Functional Modules of Olaparib Resistance Revealed by CRISPR Screens**.
It records the analysis commands and the decisions needed to reproduce candidate selection across published studies.

**Status: preparation draft, 13 October 2026.** The source notebook contains commands for six publications.
The manuscript lists ten publications. The final number of recurrence units is not yet verified.
Original MAGeCK logs, input tables, gene-summary outputs, and the original Reactome export have not been supplied to this repository.
The manuscript's reported 163 genes have not been independently reproduced here. No biological results are distributed in this draft.

## Input routes

Some studies were reanalyzed from sequencing reads; others used existing sgRNA count tables.
The analyst confirmed that some count tables were obtained directly from the original authors.
For those studies, `mageck count` was **not applicable to this reanalysis**: the analysis started with `mageck test`.
The exact authors-to-file mapping will be completed from the original records. It is not inferred from the absence of a script.

| Input available | Analysis route | Provenance to record |
| --- | --- | --- |
| FASTQ files | `mageck count`, then `mageck test` | Run accessions, FASTQ-to-sample mapping, library sequences and mapping parameters |
| Public sgRNA count table | Prepare format if necessary, then `mageck test` | Dataset/file URL, original filename, sample mapping and any transformations |
| Author-supplied sgRNA count table | Prepare format if necessary, then `mageck test` | Contributing study, original filename, date received, sample definitions and availability |

Count tables are guide-level inputs. Gene-level MAGeCK summary tables are outputs.
Existing count-table use does not mean using the original authors' selected hit lists.
Each contributing dataset must be traced to the actual input and analysis route before the final Methods are written.

## Contents

| Location | Contents |
| --- | --- |
| `scripts/mageck/` | Nine Bash scripts: four count commands and sixteen test commands, grouped by study and stage |
| `metadata/studies.tsv` | Ten publications listed in the manuscript and current input-source status |
| `metadata/samples.tsv` | Sample labels, available SRR mappings and count-table columns from the notebook |
| `metadata/comparisons.tsv` | Sixteen known comparisons and four study placeholders; all integration decisions remain pending |
| `scripts/collect_hits.py` | A new, configurable helper for positive-hit selection and recurrence counting |
| `metadata/reactome_analysis.json` | Fields to document the original enrichment analysis; unverified values remain null |
| `metadata/pathway_to_module.tsv` and `metadata/gene_to_module.tsv` | Empty schemas for the actual manual assignments |
| `provenance/command_records.json` | Source commands, paragraph references, source checksum and exploratory/final distinction |
| `docs/analysis_decisions.md` | Outstanding questions affecting reproducibility |
| `docs/reviewer_requirements.md` | Where each methodological request will be addressed |

## Recorded MAGeCK commands

The Bash scripts preserve the notebook's MAGeCK arguments and output prefixes.
Only portability wrappers, output-directory creation and existing-output checks were added.
Absolute Windows/WSL working-directory changes and environment activation were replaced by a working-directory argument.
No statistical parameters, normalization options or CNV corrections were silently added.
For Lin, the notebook's final count command was retained in the executable script; two earlier count commands remain in provenance.

Before running, place the real input files in a fresh study working directory with the relative layout used in the commands.
Activate the MAGeCK version established from the original logs. The example paths below are placeholders.

```bash
# FASTQ route: input files and the guide library must already exist.
bash scripts/mageck/lau_count.sh /absolute/path/to/1LAU
bash scripts/mageck/lau_test.sh /absolute/path/to/1LAU

# Existing count-table route: no count-generation step is needed.
bash scripts/mageck/tuano_test.sh /absolute/path/to/3TUANO
```

The recorded commands are not yet scientifically approved for a final rerun.
In particular, the Lau replicate exclusion, Fang mixed-dose test and Juhasz library handling need confirmation.
Do not treat the number of samples on a command line as proof of independent biological replication.

## Positive hits and recurrence

The postprocessing helper was written while preparing this repository; it is not asserted to be the code used in the original analysis.
Its purpose is to make the final, verified selection rule executable and auditable.
It does not force the output to contain 163 genes.

1. Replace the four study placeholders in `metadata/comparisons.tsv` with the actual comparisons.
2. Confirm input/output filenames, model identifiers, versions, normalization and replicate structure from the logs.
3. Set each comparison's `include_in_recurrence` to `yes` or `no`, with a `decision_note`.
4. For included comparisons, set `review_status` to `verified` and assign the agreed `recurrence_unit_id`.
   Multiple timepoints or sublibraries should share a unit only when that matches the explicitly documented integration rule.
5. Run the helper with the verified LFC policy. The manuscript specifies FDR < 0.05 and recurrence in at least two units;
   whether a separate LFC filter was applied still requires confirmation.

```bash
python3 scripts/collect_hits.py \
  --manifest metadata/comparisons.tsv \
  --data-root /absolute/path/containing/study-directories \
  --outdir results/verified_run \
  --fdr 0.05 --min-screens 2 \
  --lfc-policy report-only
```

`report-only` uses `pos|fdr` for selection and reports `pos|lfc` without adding a sign filter.
`require-positive` additionally requires `pos|lfc > 0`.
The choice is mandatory and must reflect the analysis being reported; the example above is not a recovered historical setting.
The helper rejects pending decisions, missing required inputs, invalid numeric fields and duplicate gene IDs within a comparison.
It uses unique recurrence units rather than counting rows, retains individual comparison statistics, and reports publication support separately.
A study with a supplied count table must not be excluded merely because a count-generation script is absent.

## Reactome and manual grouping

See `docs/reactome_workflow.md`. This draft does not rerun enrichment or invent module assignments.
The seven enrichment-based modules and literature-informed eighth group will be documented separately.

## Software and verification

The new postprocessing helper uses Python 3.9 or later and the standard library.
The historical MAGeCK version and environment are pending original logs; no unverified environment lock file is provided.

```bash
python3 -m unittest discover -s tests -v
```

The tests use small synthetic records, never study results. They check direction-specific selection,
recurrence deduplication, explicit decision requirements and failure on missing inputs.
Passing these tests does not validate the original screens or reproduce the manuscript.

## Source and citation records

Preparation sources: `Codes for others data.docx`, `ijms-4599816(1).docx`,
the supplied reviewer comments, and the analyst's clarification on author-supplied count tables.
Study DOIs are recorded in `metadata/studies.tsv`; bibliographic and sample details still need the planned provenance review.
The original manuscript, reviewer correspondence and personal correspondence are not bundled here.
The manuscript's Data Availability Statement and the release citation will be completed with actual file links and a versioned archive.
