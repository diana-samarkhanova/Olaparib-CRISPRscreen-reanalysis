#!/usr/bin/env bash
# Lau et al. 2020; mageck count.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/1LAU\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# lau_count_01
if [[ -e counts/lau_counts.count.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck count \
--norm-method median \
-l library/Human_GeCKOv2_Library_combine.csv \
--sample-label Control,OL_1,OL_2 \
--fastq \
fastq/SRR10947151.fastq.gz \
fastq/SRR10947114.fastq.gz \
fastq/SRR10947113.fastq.gz \
-n counts/lau_counts
