#!/usr/bin/env bash
# Transcribed from Codes for others data.docx; original execution logs pending.
# Study: Fang et al. 2019. Stage: mageck test.
# Paths are relative to the study working directory supplied as the first argument.
# No normalization, statistical, or CNV options were added to the source commands.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/5FANG\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# fang_test_01; source paragraphs 393-397
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

# fang_test_02; source paragraphs 400-404
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

# fang_test_03; source paragraphs 408-412
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
