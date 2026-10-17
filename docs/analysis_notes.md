# Dataset records

This page summarizes the RRA comparisons in the command records. MLE models are described in [MLE analysis](mle_analysis.md); repository coverage is listed in [data status](data_status.md).

## Recorded comparisons

| Study | Analysis in the records | Unresolved metadata |
| --- | --- | --- |
| Lau | Both OL columns counted; OL_1 tested against Control | OL_2 exclusion rationale and Control identity |
| Lin | Two samples per condition at T9; three at T12, T15 and T18. T0 is included in counting but absent from the RRA tests. | Replicate structure, DMSO3 selection and recurrence-unit definition |
| Tuano | HMLE, mesHMLE, B80-T17, B80-T5 and K5+/K19- analyzed separately | Count-table conversion, CRISPR subset, replicate structure and recurrence-unit definition |
| Juhasz | GeCKO A and B analyzed separately; the notes identify HeLa | Difference from the models/library listed in manuscript Table 1; treatment of library parts in recurrence analysis |
| Fang | Separate 5 and 10 dose comparisons and a combined-dose test | Comparison used for recurrence analysis |
| Olivieri | Olaparib_T15 versus NT_T15 | NT control identity and treatment-duration metadata |
| Zimmermann, Clements, Tsujino, Zelceski | RRA commands are not archived | Input paths and executed test commands |

Lau's recorded zero-count guide fractions are 81.90% for OL_1 and 83.07% for OL_2, with Gini indices of 0.878 and 0.888. The records do not document a QC criterion for excluding OL_2.

B80-T5 and B80-T17 are distinct cell lines in [Tuano et al.](https://doi.org/10.1186/s13059-023-02898-w). T5 and T17 are part of their names, not sampling days.

## Analysis settings

The historical RRA commands omit explicit normalization and CNV options. Their effective settings and MAGeCK version require the corresponding run logs, which are not archived. `metadata/comparisons.tsv` retains separate fields for verified normalization, replicate structure and inclusion decisions.

The repository's RRA postprocessing uses positive-selection statistics. Its application to the manuscript gene list, any additional LFC filter and the use of publication-derived hit lists remain undocumented.

## Recurrence units

The comparison manifest assigns recurrence units independently of publication IDs and individual tests. Comparisons sharing a recurrence unit count once per gene. For timepoints, doses and library parts, the aggregation rule and inclusion decision are recorded in `decision_note`.

Differences in library coverage affect which genes can be detected across studies. Gene identifiers, the tested gene sets and any union/intersection rule are therefore part of the integration record.

## References

- [MAGeCK usage](https://sourceforge.net/p/mageck/wiki/usage/)
- [MAGeCK questions and answers](https://sourceforge.net/p/mageck/wiki/QA/)
