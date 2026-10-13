# Reactome analysis and module assignment

The enrichment analysis and manual grouping are separate steps. The current repository contains their metadata schemas,
not new enrichment results or invented assignments.

## Recover the original enrichment analysis

Populate `metadata/reactome_analysis.json` from the actual analysis record:

- Exact input genes and identifier type; retain mapped and unmapped identifiers.
- Reactome release, analysis date, organism and service or software used.
- Background/universe actually used, with a file or precise definition.
- Statistical test, multiple-testing correction and FDR threshold.
- Options affecting projection, interactors and identifier mapping.
- Complete exported results before significance filtering.

If the analysis was performed on the Reactome website, record that fact and retain its exported results.
An API script may be added for a documented rerun after the settings are recovered.
A newly written script is not evidence that the historical analysis used that script.
Running against a newer Reactome release may change mappings, pathway counts and significance.

The background needs particular attention because the input studies include focused libraries.
First recover the reported analysis; document and distinguish any revised or sensitivity analysis with a different universe.
The manuscript's values (163 input genes, 733 pathways and 109 significant pathways) are presently unverified targets for comparison,
not constants that a script should force.

## Record manual grouping

In `metadata/pathway_to_module.tsv`, record the stable Reactome pathway ID, pathway name,
module ID, assignment basis, source reference and curation note.
In `metadata/gene_to_module.tsv`, record each gene-module relationship and the supporting pathway or literature reference.
Use `assignment_basis` to distinguish `reactome_enrichment` from `literature_informed` assignments.

There are seven enrichment-based modules in the manuscript framework and a literature-informed eighth group.
Do not describe this as unsupervised clustering unless such an algorithm was actually used.
Do not assume every gene lacking a significant pathway belongs to a mechanistically coherent eighth module.
Preserve the actual reasoning and source references for those assignments.

Before filling the tables, document whether multiple module memberships are allowed.
For multi-membership, retain one row per relationship and distinguish relationship counts from unique gene counts.
If the figures require one primary assignment, document the prioritization rule and retain secondary memberships in the evidence table.
The `is_primary_assignment` column is a record of a decision, not a reason to choose one automatically.

After those tables are complete, a short script can apply the curated assignments to regenerate tables and figures.
It should preserve the provenance of manual choices and check for unassigned genes, conflicting primary assignments and duplicate relationships.

Reference: https://reactome.org/userguide/analysis
