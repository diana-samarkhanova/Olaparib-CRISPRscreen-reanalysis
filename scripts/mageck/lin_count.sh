#!/usr/bin/env bash
# Lin et al. 2024; mageck count.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/2LIN\n' "$0" >&2
  exit 2
fi
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# lin_count_03
if [[ -e counts/lin_counts_OLA_DMSO3.count.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- counts
mageck count \
-l library/library_sequences_mageck.tsv \
--sample-label T0_Screen3_C,T0_Screen3_B,T0_Screen3_A,T9_OLA_B,T9_OLA_A,T9_DMSO3_B,T9_DMSO3_A,T12_OLA_C,T12_OLA_B,T12_OLA_A,T12_DMSO3_C,T12_DMSO3_B,T12_DMSO3_A,T15_OLA_C,T15_OLA_B,T15_OLA_A,T15_DMSO3_C,T15_DMSO3_B,T15_DMSO3_A,T18_OLA_C,T18_OLA_B,T18_OLA_A,T18_DMSO3_C,T18_DMSO3_B,T18_DMSO3_A \
--fastq \
fastq/SRR26388029.fastq.gz \
fastq/SRR26388030.fastq.gz \
fastq/SRR26388031.fastq.gz \
fastq/SRR26388023.fastq.gz \
fastq/SRR26388024.fastq.gz \
fastq/SRR26388027.fastq.gz \
fastq/SRR26388028.fastq.gz \
fastq/SRR26388118.fastq.gz \
fastq/SRR26388119.fastq.gz \
fastq/SRR26388120.fastq.gz \
fastq/SRR26388124.fastq.gz \
fastq/SRR26388125.fastq.gz \
fastq/SRR26388126.fastq.gz \
fastq/SRR26388104.fastq.gz \
fastq/SRR26388105.fastq.gz \
fastq/SRR26388107.fastq.gz \
fastq/SRR26388111.fastq.gz \
fastq/SRR26388112.fastq.gz \
fastq/SRR26388113.fastq.gz \
fastq/SRR26388097.fastq.gz \
fastq/SRR26388098.fastq.gz \
fastq/SRR26388099.fastq.gz \
fastq/SRR26388071.fastq.gz \
fastq/SRR26388072.fastq.gz \
fastq/SRR26388073.fastq.gz \
-n counts/lin_counts_OLA_DMSO3
