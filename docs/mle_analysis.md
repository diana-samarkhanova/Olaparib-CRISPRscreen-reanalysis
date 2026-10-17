# MAGeCK MLE analysis

The Lau, Lin, Tuano and Clements workflows use raw sgRNA count tables with `mageck mle`. This page describes their inputs, models and outputs. Run-log and result availability is listed in [data status](data_status.md).

## Environment and launch

The environment specification pins MAGeCK 0.5.9.5, Python 3.10 and NumPy 1.23.5:

```bash
conda env create -f environment-mle.yml
conda activate mageck_env
```

An existing environment with these versions can also be used. The YAML file specifies the installation dependencies; a complete export of the original environment is not archived.

Copy the required script from `scripts/mle/` into the directory containing its count table. From that directory, run:

| Study | Input file beside the script | Command |
| --- | --- | --- |
| Lau | `lau_counts.count.txt` | `python run_lau_mle.py` |
| Lin | `lin_counts_OLA_DMSO3.count.txt` | `python run_lin_mle_joint.py` |
| Tuano | `tuano_count_edited.txt` | `python run_tuano_mle.py` |
| Clements | `Clements_HeLa_WT.count.txt` and `Clements_HeLa_BRCA2KO.count.txt` | `python run_clements_mle.py` |

Adding `--dry-run` prepares and validates the inputs without fitting MLE. Each run writes a new directory under `MAGeCK_MLE_results/`, with the selected counts, design matrix, command and run metadata. The original count table is kept intact. Clements also accepts `--data-dir`; Lin and Tuano accept `--count-table` (their output still goes beside the script).

## Models and parameters

| Study | Model | Normalization |
| --- | --- | --- |
| Lau | OL_1 versus Control; one sample column in each group | total |
| Lin | All 25 samples: three T0, 11 DMSO and 11 OLA; categorical day effects plus one shared OLA-versus-DMSO effect | median |
| Tuano | Separate OL-versus-DMSO fits for HMLE, mesHMLE, B80_T5, B80_T17 and K5+/K19-; two samples per condition | median |
| Clements | Separate OL-versus-DMSO fits for HeLa WT and HeLa BRCA2KO; two samples per condition | median |

All scripts default to two threads and two permutation rounds. Clements accepts `--threads` and `--permutation-rounds`. Other options use MAGeCK defaults. The scripts do not apply CNV correction.

[The comparison manifest](../metadata/mle_comparisons.tsv) lists sample labels and settings for the nine fits in the current scripts. [Design matrices](../metadata/mle_designs/) were generated directly from those scripts. Each selected count table is written in the same sample order.

Lau uses OL_1 versus Control, matching the recorded RRA comparison, with total normalization following the recorded fallback on sparse counts. The comparison has one sample per group and no separate T0. The OL_2 exclusion rationale and the identity of Control are not documented.

Lin uses categorical day effects and one common drug coefficient across T9-T18. The model has no treatment-by-day interaction, replicate-block term or longitudinal covariance structure. The per-timepoint configuration is archived in `scripts/mle/previous/run_lin_mle.py`. For recurrence counting, multiple timepoints from the same screen contribute one unit.

Tuano's B80-T5 and B80-T17 are distinct cell lines and are fitted separately. Their recurrence-unit definition is separate from the MLE model. Tuano, Lin joint and Clements scripts use reversible aliases for whitespace or double quotes in guide/gene IDs to support the MAGeCK parser. Counts and gene groupings are preserved, and derived CSVs restore the original gene labels.

Clements inputs come from the published knockout-screen supplements: MOESM4 for HeLa WT and MOESM5 for HeLa BRCA2KO. On the worksheet `sgRNA scores and read counts`, columns A/B provide sgRNA/Gene identifiers and N/P/O/Q map to DMSO_1/DMSO_2/OL_1/OL_2. Only individual-sample read counts enter MLE. Both exports contain 77,441 guides. Reference SHA-256 hashes identify the exports, and the script validates integer-valued counts before fitting. The separate genotype fits do not test a genotype-by-drug interaction.

## Outputs

The full `*.gene_summary.txt` retains both directions of effect. The scripts also write a complete `*_ranked.csv` and a positive-hit file using:

- `treatment_vs_control|beta > 0`;
- `treatment_vs_control|fdr < 0.05` (permutation FDR, not `wald-fdr`).

`delta_beta` is the directly fitted treatment-versus-control coefficient. Clements also exports negative hits at beta < 0 and the same FDR threshold. FDR is adjusted within each fit.

`scripts/collect_hits.py` accepts RRA columns. MLE recurrence aggregation is not implemented.

## Version record

[The script record](../provenance/mle_scripts.json) lists current SHA-256 hashes, source-version hashes and validation records. The 17 October editorial update revises documentation and B80 cell-line labels; sample selection, design matrices and analysis parameters are unchanged.

Both Clements exports passed input validation and the MAGeCK 0.5.9.5 count loader with all guides and counts retained. Dry runs validate input preparation; convergence and numerical warnings require the completed fit and its log.

## Sources

- [MAGeCK input and design matrices](https://sourceforge.net/p/mageck/wiki/input/)
- [MAGeCK options](https://sourceforge.net/p/mageck/wiki/usage/)
- [MAGeCK output columns](https://sourceforge.net/p/mageck/wiki/output/)
- [Clements et al. (2020)](https://doi.org/10.1038/s41467-020-19961-w)
- [Tuano et al. (2023)](https://doi.org/10.1186/s13059-023-02898-w)

