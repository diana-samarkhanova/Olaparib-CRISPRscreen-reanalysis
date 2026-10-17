# Reactome analysis

The repository separates pathway enrichment from the subsequent manual assignment of functional modules. The original enrichment export and settings are not archived; current availability is listed in [data status](data_status.md).

## Enrichment record

`metadata/reactome_analysis.json` stores the input gene list, identifier type, organism, analysis date, Reactome release, software or service, background, statistical test, multiple-testing correction and FDR threshold. It also records projection/interactor settings and paths to the full result table and unmapped identifiers.

Each analysis is associated with its own settings and export. A rerun with another release or background is a separate analysis because pathway membership and significance can change. Library coverage is relevant to the background choice, particularly for focused screens.

## Functional modules

Pathways were grouped manually into seven modules. An eighth group was defined using the literature.

| File | Contents |
| --- | --- |
| `metadata/pathway_to_module.tsv` | Pathway ID, pathway name, module and assignment basis |
| `metadata/gene_to_module.tsv` | Gene, module, supporting pathway or reference, and primary assignment |

The `assignment_basis` values distinguish `reactome_enrichment` from `literature_informed` assignments. Each gene-module relationship occupies one row. Unique gene counts are reported separately from membership counts when a gene belongs to multiple modules. A primary assignment for a figure is recorded alongside the other memberships and its selection rationale.

Reference: [Reactome analysis guide](https://reactome.org/userguide/analysis).
