# Analysis decisions to verify

This record identifies questions that affect reproducibility. It does not select a preferred biological result.
Resolve them from the original logs, inputs and analysis records before finalizing integration.

## Data origin and scope

- Distinguish publications, experimental screens, analysis comparisons and recurrence units.
- The manuscript lists ten publications. Commands are currently available for six.
  Missing `mageck count` commands are expected when an existing count table was supplied.
  Record the exact study and original filename for author-supplied tables; do not invent a FASTQ-processing step.
- Record raw versus pre-normalized counts and every formatting or identifier transformation.
  `tuano_count_edited.txt` requires its original table and the changes used to create it.
- Confirm which knockout datasets were selected from studies containing KO, CRISPRi or CRISPRa experiments.
- Preserve focused-library descriptions and their target universes. Not all source datasets are genome-wide.

## MAGeCK settings

The source commands use `mageck test` and do not explicitly set normalization or CNV correction.
Record the installed version, effective defaults, filtering, replicate handling, normalization,
statistical adjustment and any options shown in the actual logs. Do not infer that omitted options meant no normalization.
Do not substitute a version remembered from another project.

Separate the final hit-selection FDR threshold from any internal MAGeCK algorithm parameters.
The manuscript says FDR < 0.05. Verify that the positive-selection field was used,
whether an additional LFC criterion was used, and whether any original-publication hit lists entered the final set.
Using data from a sensitivity-oriented publication does not itself mean using negative-selection hits.

## Study-specific questions

| Study | Evidence in the supplied notebook | Decision or record needed |
| --- | --- | --- |
| Lau | Both OL replicates counted; only OL_1 tested against one control | Establish the reason for OL_2 exclusion and assess its consequence. Report the actual single-treatment/single-control design if retained. |
| Lin | T9 has two treatment/control samples; T12, T15 and T18 have three each. T0 is counted but not used in the shown tests. | Verify biological replicate definitions, DMSO3 control choice, and how timepoints contribute to one recurrence unit. |
| Tuano | Five tests using HMLE, mesHMLE, B80_T17, B80_T5 and K5+/K19- labels | Verify cell-model identities and the meaning of T5/T17; do not assume these labels denote treatment timepoints. Record count-table edits and the KO subset. |
| Juhasz | Separate GeCKO A and B count/test commands. Notebook text says HeLa; manuscript Table 1 lists different models and library. | Reconcile models and library using source accessions; decide whether/how A/B were combined and prevent duplicate recurrence support. |
| Fang | Separate Olap_5 and Olap_10 tests and one combined-dose test | Identify which outputs entered integration. Distinguish dose groups from identical-condition biological replicates. |
| Olivieri | Olaparib_T15 versus NT_T15 | Confirm timepoint, NT control and source-table provenance; reconcile with the manuscript table. |
| Zimmermann, Clements, Tsujino, Zelceski | Listed in manuscript; their commands are not in the current notebook | Add the actual supplied/public count tables and test logs or commands. Lack of FASTQ code is not a reason for exclusion. |

For Lau, the recorded zero-guide proportions are 81.90% for OL_1 and 83.07% for OL_2,
with Gini indices 0.878 and 0.888. These summaries alone do not establish why only OL_2 should be excluded.
The count of mapped reads and a Gini index are QC observations, not a pre-established exclusion rule.

## Recurrence and enrichment

- Document union/intersection or any other method used to combine results before recurrence counting.
- Define the exact recurrence unit. Sharing a publication or model is not by itself a complete independence definition.
- Keep each comparison's FDR/LFC and source labels; do not assign one pooled FDR to a recurrent gene without a defined meta-analysis.
- Verify identifier conversions, duplicate resolution and eligibility when libraries cover different genes.
- Recover the original Reactome release, background, input list and complete export.
- Keep manual pathway grouping distinct from statistical enrichment; record overlap handling explicitly.

## Release readiness

The final code should reproduce the revised, verified results. If correcting a method changes the gene count,
update the manuscript and supplementary material consistently instead of tuning code to restore 163.
Complete the source links, availability descriptions, software versions and release citation after the records are checked.

Official technical references consulted while preparing this draft:

- https://sourceforge.net/p/mageck/wiki/usage/
- https://sourceforge.net/p/mageck/wiki/QA/
- https://reactome.org/userguide/analysis
