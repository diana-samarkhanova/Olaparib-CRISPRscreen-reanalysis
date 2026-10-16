# Olaparib CRISPR screen reanalysis

Code and dataset records for *Beyond Homologous Recombination: Functional Modules of Olaparib Resistance Revealed by CRISPR Screens*.

The RRA commands cover six of the ten publications in the review. MLE scripts are available for Lau, Lin, Tuano and Clements. Original run logs, count tables and analysis outputs are being added. The reported set of 163 recurrent genes has not yet been reproduced from this repository.

## Data sources

FASTQ datasets were processed with `mageck count`, followed by `mageck test`. Public count tables and tables obtained directly from the authors entered the analysis at `mageck test`; no count-generation step was needed for those datasets.

Count tables contain sgRNA read counts, not preselected gene lists. Sources, filenames and sample labels are recorded in `metadata/`. The remaining author-supplied files will be identified when the run records are added.

## Files

| Path | Contents |
| --- | --- |
| `scripts/mageck/` | RRA count and test commands, grouped by study |
| `scripts/mle/` | MLE scripts for Lau, Lin, Tuano and Clements |
| `metadata/mle_comparisons.tsv` | MLE samples, settings and design-matrix paths |
| `metadata/mle_designs/` | Design matrices generated from the current MLE scripts |
| `scripts/collect_hits.py` | Positive-hit selection and recurrence counting |
| `metadata/` | Studies, sample labels, comparisons and Reactome settings |
| `provenance/command_records.json` | Original commands and their locations in the analysis notes |
| `docs/` | Dataset questions and Reactome methods |
| `tests/` | Tests for hit selection and recurrence counting |

The shell scripts in `scripts/mageck/` preserve the original arguments and output prefixes. Each takes the study working directory as its only argument. Lin's final count command is used; the two earlier versions remain in the command record.

## Running RRA

Activate the MAGeCK environment and place the input files in the relative paths shown in each script. The original software version and effective defaults still need to be checked against the logs.

```bash
# From FASTQ files
bash scripts/mageck/lau_count.sh /path/to/1LAU
bash scripts/mageck/lau_test.sh /path/to/1LAU

# From an existing count table
bash scripts/mageck/tuano_test.sh /path/to/3TUANO
```

Use a fresh output directory for each run. See [analysis notes](docs/analysis_notes.md) for the comparisons that need clarification.

## Running MLE

See [MLE analysis](docs/mle_analysis.md) for the environment, inputs, commands and model definitions. Lin uses one shared OLA effect adjusted for sampling day. Tuano cell models and Clements genotypes are fitted separately. The earlier Lin per-timepoint script is retained in `scripts/mle/previous/`.

MLE positive-hit files use `treatment_vs_control|beta > 0` and permutation `treatment_vs_control|fdr < 0.05`. The recurrence script below currently accepts RRA outputs; it does not combine MLE results.

## Selecting recurrent genes from RRA

`collect_hits.py` is a new postprocessing script. Its selection rules will be checked against the original gene lists.

Before running it, complete `metadata/comparisons.tsv`: replace the missing-study rows, record inclusion decisions, and assign a `recurrence_unit_id` to each included comparison. Included rows also require a model ID, `review_status=verified`, and a decision note. Comparisons assigned to the same recurrence unit contribute only once per gene.

Example invocation:

```bash
python3 scripts/collect_hits.py \
  --manifest metadata/comparisons.tsv \
  --data-root /path/to/study-directories \
  --outdir results/run_01 \
  --fdr 0.05 --min-screens 2 \
  --lfc-policy report-only
```

`report-only` selects by `pos|fdr` and retains `pos|lfc` in the output. `require-positive` also requires `pos|lfc > 0`. The LFC rule used in the original analysis remains to be confirmed; the example above is not a final analysis setting.

Outputs include per-comparison evidence, hit counts, recurrent genes and their supporting studies. The script requires Python 3.9 or later and uses only the standard library.

## Reactome

The enrichment settings and module-assignment tables are awaiting the original export. Seven modules were based on enrichment results; the eighth group was informed by the literature. See [Reactome methods](docs/reactome_workflow.md) for the information to record.

## Tests

```bash
python3 -m unittest discover -s tests -v
```

Tests use synthetic data to check positive-hit selection, recurrence counting and input validation. Study datasets have not yet been rerun here.

