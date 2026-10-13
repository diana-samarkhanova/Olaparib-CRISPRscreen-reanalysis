#!/usr/bin/env bash
# Fang et al. 2019; mageck test.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/5FANG\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# fang_test_01
if [[ -e results/Olap5_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- results
mageck test \
-k count/all.count.txt \
-t Olap_5.1,Olap_5.2 \
-c DMSO.1,DMSO.2 \
-n results/Olap5_vs_DMSO

# fang_test_02
if [[ -e results/Olap10_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- results
mageck test \
-k count/all.count.txt \
-t Olap_10.1,Olap_10.2 \
-c DMSO.1,DMSO.2 \
-n results/Olap10_vs_DMSO

# fang_test_03
if [[ -e results/Olap_all_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- results
mageck test \
-k count/all.count.txt \
-t Olap_5.1,Olap_5.2,Olap_10.1,Olap_10.2 \
-c DMSO.1,DMSO.2 \
-n results/Olap_all_vs_DMSO
