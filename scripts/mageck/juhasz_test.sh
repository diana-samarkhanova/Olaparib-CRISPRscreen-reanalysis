#!/usr/bin/env bash
# Juhasz et al. 2020; mageck test.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/4JUHASZ\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# juhasz_test_01
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

# juhasz_test_02
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
