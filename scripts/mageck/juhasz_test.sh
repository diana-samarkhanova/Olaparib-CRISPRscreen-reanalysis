#!/usr/bin/env bash
# Transcribed from Codes for others data.docx; original execution logs pending.
# Study: Juhasz et al. 2020. Stage: mageck test.
# Paths are relative to the study working directory supplied as the first argument.
# No normalization, statistical, or CNV options were added to the source commands.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/4JUHASZ\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# juhasz_test_01; source paragraphs 330-334
if [[ -e counts/GeCKO_A_OL_vs_DMSO_test.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck test \
-k counts/GeCKO_A_OL_VS_DMSO.count.txt \
-t A_OL_1,A_OL_2 \
-c A_DMSO_1,A_DMSO_2 \
-n counts/GeCKO_A_OL_vs_DMSO_test

# juhasz_test_02; source paragraphs 361-365
if [[ -e counts/GeCKO_B_OL_vs_DMSO_test.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck test \
-k counts/GeCKO_B_OL_VS_DMSO.count.txt \
-t B_OL_1,B_OL_2 \
-c B_DMSO_1,B_DMSO_2 \
-n counts/GeCKO_B_OL_vs_DMSO_test
