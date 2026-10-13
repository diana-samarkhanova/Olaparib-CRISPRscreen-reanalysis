#!/usr/bin/env bash
# Transcribed from Codes for others data.docx; original execution logs pending.
# Study: Lau et al. 2020. Stage: mageck test.
# Paths are relative to the study working directory supplied as the first argument.
# No normalization, statistical, or CNV options were added to the source commands.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/1LAU\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# lau_test_01; source paragraphs 45-49
if [[ -e counts/lau_test_Control_vs_OL1.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck test \
-k counts/lau_counts.count.txt \
-t OL_1 \
-c Control \
-n counts/lau_test_Control_vs_OL1
