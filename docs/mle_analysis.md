# MAGeCK MLE analysis

MLE scripts are provided for four studies. They use sgRNA count tables and call `mageck mle`; the RRA scripts use `mageck test`. Full MLE run logs and result tables have not yet been archived here. This addition documents the code and settings, rather than a completed reanalysis of every study in the review.

## Environment and launch

The environment specification pins MAGeCK 0.5.9.5, Python 3.10 and NumPy 1.23.5:

```bash
conda env create -f environment-mle.yml
conda activate mageck_env
```

Use the existing environment if these versions are already installed. The YAML file is an installation specification, not a complete export of the original environment.

Copy the required script from `scripts/mle/` into the directory containing its count table. From that directory, run:

| Study | Input file beside the script | Command |
| --- | --- | --- |
| Lau | `lau_counts.count.txt` | `python run_lau_mle.py` |
| Lin | `lin_counts_OLA_DMSO3.count.txt` | `python run_lin_mle_joint.py` |
| Tuano | `tuano_count_edited.txt` | `python run_tuano_mle.py` |
| Clements | `Clements_HeLa_WT.count.txt` and `Clements_HeLa_BRCA2KO.count.txt` | `python run_clements_mle.py` |

Adding `--dry-run` prepares and validates the inputs without fitting MLE. Each run writes a new directory under `MAGeCK_MLE_results/`, with the selected counts, design matrix, command and run metadata. The original count table is kept intact. Clements also accepts `--data-dir`; Lin and Tuano accept `--count-table` (their output still goes beside the script).

## Models and parameters

| Study | Model | Requested normalization |
| --- | --- | --- |
| Lau | OL_1 versus Control; one sample column in each group | median; sparse inputs can block the run |
| Lin | All 25 samples: three T0, 11 DMSO and 11 OLA; categorical day effects plus one shared OLA-versus-DMSO effect | median |
| Tuano | Separate OL-versus-DMSO fits for HMLE, mesHMLE, B80_T5, B80_T17 and K5+/K19-; two samples per condition | median |
| Clements | Separate OL-versus-DMSO fits for HeLa WT and HeLa BRCA2KO; two samples per condition | median |

All scripts default to two threads and two permutation rounds. Clements allows these to be changed with `--threads` and `--permutation-rounds`. Other MAGeCK options use the defaults of the installed version; no CNV correction is supplied.

[The comparison manifest](../metadata/mle_comparisons.tsv) lists sample labels and settings for the nine fits in the current scripts. [Design matrices](../metadata/mle_designs/) were generated directly from those scripts. Each selected count table is written in the same sample order.

Lau's sample selection follows the recorded RRA comparison. OL_2 is not selected; the reason for this exclusion remains unresolved. The earlier script requested total following a reported median-to-total fallback. The 17 October 2026 version requests median and refuses inputs that would trigger the known sparsity fallback; a successful new run is required before its results can be described as median-normalized. The label `Control` is retained without assuming that it denotes DMSO. This fit has no within-group replication and no separate T0 sample.

All four current scripts check selected samples before fitting, including during `--dry-run`, and reject any sample with more than 45% zero-count guides. This is the automatic total-fallback threshold in MAGeCK 0.5.9.5. After each fit they check the engine log before exporting ranked or hit CSVs, and save `*.normalization.json`. A warning at more than 30% but at most 45% zeros is not itself a fallback. The `normalization` field in `run_metadata.json` and the manifest records the requested method; it is not proof that a run completed. See [normalization policy](normalization.md).

Lin's shared drug coefficient assumes no treatment-by-day interaction. It is not a slope over time. A/B/C sample labels are not used to infer pairing or correlation between repeated cultures. The earlier separate-timepoint script is retained as `scripts/mle/previous/run_lin_mle.py` to reproduce that alternative configuration; it is not the current joint model. Timepoints from one screen must not be counted as independent recurrence evidence.

Tuano's B80-T5 and B80-T17 are distinct cell lines, not sampling timepoints; they are fitted separately. The recurrence unit must be defined explicitly before integration. This corrects the earlier metadata labels and does not change selected samples or design matrices. The current Tuano script includes the identifier-parser fix: whitespace or double quotes in sgRNA/gene IDs are replaced by reversible, collision-free aliases for MAGeCK. The same handling is present in the Lin joint and Clements scripts. These aliases do not change counts or merge genes; derived CSVs restore original gene labels.

Clements inputs were exported from the published knockout-screen supplements: MOESM4 for HeLa WT and MOESM5 for HeLa BRCA2KO. On the worksheet `sgRNA scores and read counts`, columns A/B provide sgRNA/Gene identifiers; N/P/O/Q become DMSO_1/DMSO_2/OL_1/OL_2. No supplied RRA statistic columns were used. Both saved exports contain 77,441 guides. The script records their reference SHA-256 hashes and validates integer-valued counts before fitting. The separate genotype fits do not test a genotype-by-drug interaction.

## Outputs

The full `*.gene_summary.txt` retains both directions of effect. The scripts also write a complete `*_ranked.csv` and a positive-hit file using:

- `treatment_vs_control|beta > 0`;
- `treatment_vs_control|fdr < 0.05` (permutation FDR, not `wald-fdr`).

`delta_beta` is the directly fitted treatment-versus-control coefficient. It is not a subtraction of two independently fitted effects. Clements also exports negative hits at beta < 0 and the same FDR threshold. FDR adjustment is within each fit, without adjustment across studies or comparisons.

No MLE recurrence list is produced here. `scripts/collect_hits.py` currently expects RRA columns. Inclusion in the cross-screen analysis requires the corresponding completed outputs and the study's recurrence-unit decision.

## Version record

[The script record](../provenance/mle_scripts.json) contains current SHA-256 hashes and retains the earlier hashes and commit reference. The 17 October 2026 update explicitly requests median in Lau and adds normalization checks to all current MLE scripts. Publishing this update does not rerun MAGeCK or update old result tables. The archived Lin per-timepoint script is unchanged.

The Clements exports passed input validation and the original MAGeCK 0.5.9.5 count loader without losing guides or changing counts. A dry run checks input preparation, not model convergence. Numerical warnings observed during an MLE permutation stage still require review against the completed run and result table; the screenshot alone does not identify a particular study.

## Sources

- [MAGeCK input and design matrices](https://sourceforge.net/p/mageck/wiki/input/)
- [MAGeCK options](https://sourceforge.net/p/mageck/wiki/usage/)
- [MAGeCK output columns](https://sourceforge.net/p/mageck/wiki/output/)
- [Clements et al. (2020)](https://doi.org/10.1038/s41467-020-19961-w)
- [Tuano et al. (2023), cell-line identities](https://doi.org/10.1186/s13059-023-02898-w)

