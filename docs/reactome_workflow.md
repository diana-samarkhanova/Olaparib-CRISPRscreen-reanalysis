# Reactome analysis

The original enrichment export and settings are needed to complete this part of the repository.

## Enrichment

Record the following in `metadata/reactome_analysis.json`:

- Input gene file, identifier type and mapped/unmapped identifiers.
- Reactome release, analysis date, organism and analysis service or software.
- Background definition and gene list, where applicable.
- Statistical test, multiple-testing correction and FDR threshold.
- Projection and interactor settings.
- Full result table before significance filtering.

If the original analysis used the Reactome website, retain that export and describe the web analysis. Any later scripted rerun should record its own settings and database release. Changes in release or background may change pathway membership and significance. The background is especially relevant because the source studies include focused libraries.

## Functional modules

Pathways were grouped manually into seven modules. An eighth group was defined using the literature.

| File | Required information |
| --- | --- |
| `metadata/pathway_to_module.tsv` | Reactome pathway ID, pathway name, module and basis for assignment |
| `metadata/gene_to_module.tsv` | Gene, module, supporting pathway or reference, and any primary assignment |

Use `reactome_enrichment` or `literature_informed` in `assignment_basis` to distinguish the two sources of evidence. Record the rule for genes assigned to more than one module. If multiple memberships are retained, use one row per gene-module relationship and report unique gene counts separately. If a figure uses one primary assignment, document how it was chosen and retain the other assignments in the table.

Once the assignments are complete, they can be applied by script to regenerate the tables and figures.

Reference: [Reactome analysis guide](https://reactome.org/userguide/analysis).
