# Reproducibility requirements

This table maps the supplied reviewer requests to repository records and manuscript revisions.
It records work to complete and does not claim that the revision has been performed.

| Request | Repository evidence | Manuscript location | Current status |
| --- | --- | --- | --- |
| Search date, selection and dataset accessions | `metadata/studies.tsv`, `metadata/samples.tsv`; add search record and exact source links | Methods 4.1 and supplementary source table | Partial |
| Distinguish FASTQ and supplied-count routes | Study input provenance and stage-specific scripts | Methods 4.1-4.2 | Routes documented; author-to-file mapping pending |
| MAGeCK version and exact commands/settings | `scripts/mageck/`, original logs and environment record | Methods 4.2 | Commands for six publications; logs pending |
| Replicates, normalization and treatment-control comparisons | `metadata/comparisons.tsv` and `metadata/samples.tsv` | Methods 4.2 and comparison table | Sample labels recorded; effective settings pending |
| Positive-selection threshold and direction | `scripts/collect_hits.py`, per-comparison gene summaries | Methods 4.2; reconcile Discussion 3.3 | Historical filtering decision pending |
| Precise definition of recurrence | Verified recurrence-unit IDs and explicit deduplication | Methods 4.2-4.3 and Figure 1 | Pending |
| Identify studies/models supporting each recurrent hit | Per-comparison evidence and recurrent-gene support output | Supplementary gene table | Awaiting actual results |
| Count comparisons/publications correctly | Separate comparison IDs, study IDs and recurrence units | Methods, Results and figure legends | Catalog started |
| Uniform reanalysis versus borrowed hit lists | Input type, original file and analysis provenance per dataset | Methods 4.2 | Existing counts are distinguished from selected gene lists; remaining provenance pending |
| Reactome version/background and complete results | `metadata/reactome_analysis.json` and original export | Methods 4.3 | Awaiting original files |
| Gene/module assignments and overlaps | Two assignment tables and a documented overlap rule | Methods 4.3 and supplementary assignments | Schemas prepared; assignments pending |
| Seven enriched modules and literature-informed eighth group | Assignment basis and literature sources | Results, Figure 1 and Discussion | Explicit distinction documented |
| Data/code availability | Exact input sources, final repository URL and versioned archive | Data Availability Statement | Final links pending |

The scientific claims, reference corrections and figure-label changes requested by reviewers remain separate manuscript tasks.
This project addresses the computational and provenance portion of the revision.
