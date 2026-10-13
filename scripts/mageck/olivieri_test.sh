#!/usr/bin/env bash
# Olivieri et al. 2020; mageck test.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/7OLIVIERI\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# olivieri_test_01
if [[ -e counts/dataset2_OL_vs_NT_T15_test.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck test \
-k counts/Dataset_S2_readcounts.txt \
-t S02_Olaparib_T15_A,S02_Olaparib_T15_B \
-c S02_NT_T15_A,S02_NT_T15_B \
-n counts/dataset2_OL_vs_NT_T15_test
