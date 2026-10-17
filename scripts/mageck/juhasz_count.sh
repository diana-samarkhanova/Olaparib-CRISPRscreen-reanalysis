#!/usr/bin/env bash
# Juhasz et al. 2020; mageck count.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/4JUHASZ\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# juhasz_count_01
if [[ -e counts/GeCKO_A_OL_VS_DMSO.count.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck count \
--norm-method median \
-l library/Human_GeCKOv2_Library_A_3_mageck.csv \
--sample-label A_OL_1,A_OL_2,A_DMSO_1,A_DMSO_2 \
--fastq \
fastq/SRR11144451.fastq.gz \
fastq/SRR11144452.fastq.gz \
fastq/SRR11144453.fastq.gz \
fastq/SRR11144450.fastq.gz \
-n counts/GeCKO_A_OL_VS_DMSO

# juhasz_count_02
if [[ -e counts/GeCKO_B_OL_VS_DMSO.count.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck count \
--norm-method median \
-l library/Human_GeCKOv2_Library_B_1_mageck.csv \
--sample-label B_OL_1,B_OL_2,B_DMSO_1,B_DMSO_2 \
--fastq \
fastq/SRR11144455.fastq.gz \
fastq/SRR11144456.fastq.gz \
fastq/SRR11144457.fastq.gz \
fastq/SRR11144458.fastq.gz \
-n counts/GeCKO_B_OL_VS_DMSO
