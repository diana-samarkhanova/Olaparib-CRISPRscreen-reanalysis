# Starter validation

Validation date: 13 October 2026.

- Nine study/stage Bash scripts passed `bash -n`.
- All twenty retained MAGeCK commands matched their source-notebook text exactly.
- Two earlier Lin count commands were retained in provenance, not in execution scripts.
- Five synthetic postprocessing tests passed: unique recurrence units, positive-selection fields, explicit LFC policy, pending decisions and missing/duplicated source files.
- TSV and JSON metadata passed structural checks.

MAGeCK was not executed on study data. No original count tables, gene summaries or Reactome exports were available.
These checks validate the starter code and transcription, not the manuscript gene count or scientific conclusions.
