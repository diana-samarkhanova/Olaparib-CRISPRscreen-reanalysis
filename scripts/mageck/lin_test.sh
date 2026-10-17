#!/usr/bin/env bash
# Lin et al. 2024; mageck test.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/2LIN\n' "$0" >&2
  exit 2
fi
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/median_normalization.sh"
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# lin_test_01
if [[ -e test/T9_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
--norm-method median \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T9_OLA_B,T9_OLA_A \
-c T9_DMSO3_B,T9_DMSO3_A \
-n test/T9_OLA_vs_DMSO3
require_median_log "test/T9_OLA_vs_DMSO3.log"

# lin_test_02
if [[ -e test/T12_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
--norm-method median \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T12_OLA_C,T12_OLA_B,T12_OLA_A \
-c T12_DMSO3_C,T12_DMSO3_B,T12_DMSO3_A \
-n test/T12_OLA_vs_DMSO3
require_median_log "test/T12_OLA_vs_DMSO3.log"

# lin_test_03
if [[ -e test/T15_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
--norm-method median \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T15_OLA_C,T15_OLA_B,T15_OLA_A \
-c T15_DMSO3_C,T15_DMSO3_B,T15_DMSO3_A \
-n test/T15_OLA_vs_DMSO3
require_median_log "test/T15_OLA_vs_DMSO3.log"

# lin_test_04
if [[ -e test/T18_OLA_vs_DMSO3.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mkdir -p -- test
mageck test \
--norm-method median \
-k counts/lin_counts_OLA_DMSO3.count.txt \
-t T18_OLA_C,T18_OLA_B,T18_OLA_A \
-c T18_DMSO3_C,T18_DMSO3_B,T18_DMSO3_A \
-n test/T18_OLA_vs_DMSO3
require_median_log "test/T18_OLA_vs_DMSO3.log"
