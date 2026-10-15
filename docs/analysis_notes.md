# Analysis notes

The following RRA records need checking against the original files and run logs before the metadata are finalized. MLE sample selections, normalization and design matrices are documented separately in [MLE analysis](mle_analysis.md).

## Inputs and settings

Record the source and original filename of each count table, including tables received from authors. Keep any formatting or gene-identifier changes with the original file. For studies using more than one CRISPR approach, identify the knockout subset used here.

The recorded RRA commands do not explicitly set normalization or CNV correction. Recover the MAGeCK version and effective defaults from the logs, along with replicate handling, filtering and statistical settings. The final hit-selection threshold is separate from MAGeCK's internal parameters.

Confirm that positive-selection results were used throughout, whether an additional LFC filter was applied, and whether any hit lists were taken directly from the original publications. Data from a sensitivity-focused study can still be reanalyzed for positive selection.

## Dataset details

| Study | Recorded analysis | To confirm |
| --- | --- | --- |
| Lau | Both OL replicates counted; only OL_1 tested against one control | Reason for excluding OL_2 and its effect on the result |
| Lin | Two treatment/control samples at T9; three at T12, T15 and T18. T0 counted but absent from the shown tests. | Replicate structure, DMSO3 control choice and grouping of timepoints |
| Tuano | Five comparisons: HMLE, mesHMLE, B80_T17, B80_T5 and K5+/K19- | Model identities, meaning of T5/T17, KO subset and edits to `tuano_count_edited.txt` |
| Juhasz | GeCKO A and B analyzed separately; notes describe HeLa | Model/library discrepancy with manuscript Table 1 and handling of A/B results |
| Fang | Separate 5 and 10 dose comparisons plus a combined-dose test | Which results entered integration; distinction between dose groups and replicates |
| Olivieri | Olaparib_T15 versus NT_T15 | Timepoint, control identity and source-table metadata |
| Zimmermann, Clements, Tsujino, Zelceski | Commands not yet included | Actual input tables and `mageck test` commands from the logs |

Lau's recorded zero-guide proportions are 81.90% for OL_1 and 83.07% for OL_2; Gini indices are 0.878 and 0.888. These values alone do not explain exclusion of OL_2.

## Combining results

Define the recurrence unit separately from the publication and individual comparison. Record how timepoints, doses and library parts were combined, including whether the rule used a union or intersection of hit lists. Retain each comparison's FDR and LFC.

Check gene identifiers, duplicates and differences in library coverage. Focused libraries should be identified as such. For Reactome, recover the original background, release and complete output before assessing the enrichment results.

The source record also needs the search date, inclusion/exclusion decisions, dataset accessions and links to the files used.

## References

- [MAGeCK usage](https://sourceforge.net/p/mageck/wiki/usage/)
- [MAGeCK questions and answers](https://sourceforge.net/p/mageck/wiki/QA/)
- [Reactome analysis guide](https://reactome.org/userguide/analysis)

