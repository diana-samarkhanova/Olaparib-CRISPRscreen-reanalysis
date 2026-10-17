# Normalization configuration, 17 October 2026

All current `mageck test` and `mageck mle` commands explicitly request `--norm-method median`. Count commands also specify median for their normalization/QC outputs; the raw `*.count.txt` remains the input to downstream analysis. Never substitute normalized counts for the raw input.

## Requested method and actual behavior

The target version is MAGeCK 0.5.9.5, pinned in `environment-mle.yml`. Its `mageck/mageckCountNorm.py::normalizeCounts` can replace median with total for the entire comparison if any sample has a zero median factor or more than 45% zero-count guides. A warning above 30% and up to 45% zeros only recommends reviewing normalization and does not by itself switch methods. Median normalization remains a scientific choice that must be suitable for the screen; these checks only prevent a silent method mismatch.

Current MLE scripts stop before fitting when a selected sample has more than 45% zero-count guides. This check also runs with `--dry-run`. They preserve all input guides and counts. After fitting, they require final normalization factors in the engine log and reject the explicit switch-to-total warning before exporting ranked/hit CSVs. Each checked fit writes `*.normalization.json`; read it together with the engine log and exit status.

RRA scripts use `median_normalization.sh` to check each completed test's log and return a nonzero exit status if MAGeCK switched to total or the log cannot establish normalization. Raw MAGeCK outputs may already exist after a rejected test; retain them for diagnosis but exclude that run from median-normalized analyses. Keep this helper beside the test scripts. Count-stage logs may still show a fallback used for QC; it does not change the raw count table or establish the normalization used by a later test.

The `normalization` field in MLE metadata describes the requested configuration. A dry run, changed source code, or the presence of `--norm-method median` is not evidence of a completed median-normalized analysis. The checks target MAGeCK 0.5.9.5's log messages; validate them again before upgrading the engine.

## What needs a new run

| Analysis | Change | Required evidence |
| --- | --- | --- |
| Lau MLE | `total` changed to `median`; fallback rejected | A new accepted fit and its logs, then regenerate any downstream results derived from it |
| Lin, Tuano, Clements MLE | Median was already requested; preflight and log checks added | Inspect completed-run logs, or rerun if actual normalization cannot be established |
| Published RRA scripts | Median is now explicit; fallback rejected | Inspect historical logs before attributing median to old outputs; rerun if needed |

Lau's recorded OL_1 zero-guide fraction is 81.90%. If that is also true of the current input, the new MLE script will stop. It is not valid to force the manuscript wording by dropping zero-count guides, inventing counts, or overriding MAGeCK's safeguard. Review the table's provenance and QC, then make and document an appropriate normalization decision. Total normalization is still normalization; it is a different method that must be reported accurately if used.

The original commands in `provenance/command_records.json`, the archived Lin script and existing run outputs are historical records and are not rewritten. No study was rerun for this code update. Lau's `Control` identity and OL_2 exclusion, Lin's joint day-adjusted model, MLE recurrence integration, and missing Reactome provenance still need their own resolution; this normalization change does not resolve them.
