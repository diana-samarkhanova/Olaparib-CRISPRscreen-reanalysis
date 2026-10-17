# Data availability and verification status

Status of the files in this repository as of 17 October 2026. A missing archive record does not establish whether an analysis was performed elsewhere.

| Component | Repository status |
| --- | --- |
| RRA commands | Available for Lau, Lin, Tuano, Juhasz, Fang and Olivieri |
| MLE scripts | Available for Lau, Lin, Tuano and Clements; nine current comparisons |
| Input count tables | Not included; sources and paths are recorded in `metadata/` |
| Complete run logs and gene-summary outputs | Not archived |
| RRA recurrence analysis | Script available; comparison inclusion and recurrence-unit decisions remain unresolved |
| MLE recurrence analysis | Not implemented in `collect_hits.py` |
| Reported 163 recurrent genes | Not reproduced from the archived repository files |
| Reactome analysis | Original export, release, background and analysis settings are not archived |
| Pathway and gene module assignments | Table headers available; assignment records are not populated |

The RRA manifest uses `PENDING` and `awaiting_logs` for unresolved records. These values prevent the recurrence script from silently treating unverified comparisons as included results. Missing values in `metadata/reactome_analysis.json` likewise remain `null`; the manuscript-reported counts are marked as unverified.

Dataset-specific gaps are listed in [dataset records](analysis_notes.md). Lau's control identity and OL_2 exclusion, replicate structure, effective RRA normalization, the historical LFC selection rule and cross-screen aggregation decisions require supporting records. The provenance directory preserves the original commands and file checksums.

Validation so far covers synthetic tests, MLE design matrices and input preparation. Both Clements exports, with 77,441 guides each, passed the MAGeCK 0.5.9.5 count loader without changes to counts or guide membership. Full study fits and model convergence have not been verified by these checks.
