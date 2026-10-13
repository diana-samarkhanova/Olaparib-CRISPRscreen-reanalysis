#!/usr/bin/env bash
# Transcribed from Codes for others data.docx; original execution logs pending.
# Study: Lin et al. 2024. Stage: mageck test.
# Paths are relative to the study working directory supplied as the first argument.
# No normalization, statistical, or CNV options were added to the source commands.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/2LIN\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# lin_test_01; source paragraphs 225-229
if [[ -e test/T9_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T9_OLA_B,T9_OLA_A \
-c T9_DMSO3_B,T9_DMSO3_A \
-n test/T9_OLA_vs_DMSO3

# lin_test_02; source paragraphs 231-235
if [[ -e test/T12_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T12_OLA_C,T12_OLA_B,T12_OLA_A \
-c T12_DMSO3_C,T12_DMSO3_B,T12_DMSO3_A \
-n test/T12_OLA_vs_DMSO3

# lin_test_03; source paragraphs 237-241
if [[ -e test/T15_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T15_OLA_C,T15_OLA_B,T15_OLA_A \
-c T15_DMSO3_C,T15_DMSO3_B,T15_DMSO3_A \
-n test/T15_OLA_vs_DMSO3

# lin_test_04; source paragraphs 243-247
if [[ -e test/T18_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T18_OLA_C,T18_OLA_B,T18_OLA_A \
-c T18_DMSO3_C,T18_DMSO3_B,T18_DMSO3_A \
-n test/T18_OLA_vs_DMSO3
