#!/usr/bin/env bash
# Tuano et al. 2023; mageck test.
# Usage: pass the study working directory as the first argument.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  printf 'Usage: bash %s /absolute/path/to/3TUANO\n' "$0" >&2
  exit 2
fi
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/median_normalization.sh"
cd -- "$1"
command -v mageck >/dev/null || { printf 'Activate the verified MAGeCK environment first.\n' >&2; exit 127; }

# tuano_test_01
if [[ -e HMLE_Olap_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mageck test \
--norm-method median \
-k tuano_count_edited.txt \
-t HMLE_OL_1,HMLE_OL_2 \
-c HMLE_DMSO_1,HMLE_DMSO_2 \
-n HMLE_Olap_vs_DMSO
require_median_log "HMLE_Olap_vs_DMSO.log"

# tuano_test_02
if [[ -e mesHMLE_Olap_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mageck test \
--norm-method median \
-k tuano_count_edited.txt \
-t mesHMLE_OL_1,mesHMLE_OL_2 \
-c mesHMLE_DMSO_1,mesHMLE_DMSO_2 \
-n mesHMLE_Olap_vs_DMSO
require_median_log "mesHMLE_Olap_vs_DMSO.log"

# tuano_test_03
if [[ -e B80_T17_Olap_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mageck test \
--norm-method median \
-k tuano_count_edited.txt \
-t B80_T17_OL_1,B80_T17_OL_2 \
-c B80_T17_DMSO_1,B80_T17_DMSO_2 \
-n B80_T17_Olap_vs_DMSO
require_median_log "B80_T17_Olap_vs_DMSO.log"

# tuano_test_04
if [[ -e B80_T5_Olap_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mageck test \
--norm-method median \
-k tuano_count_edited.txt \
-t B80_T5_OL_1,B80_T5_OL_2 \
-c B80_T5_DMSO_1,B80_T5_DMSO_2 \
-n B80_T5_Olap_vs_DMSO
require_median_log "B80_T5_Olap_vs_DMSO.log"

# tuano_test_05
if [[ -e K5_K19_Olap_vs_DMSO.gene_summary.txt ]]; then
  printf 'Output already exists; use a new study working directory to retain provenance.\n' >&2
  exit 2
fi
mageck test \
--norm-method median \
-k tuano_count_edited.txt \
-t K5+/K19-_OL_1,K5+/K19-_OL_2 \
-c K5+/K19-_DMSO_1,K5+/K19-_DMSO_2 \
-n K5_K19_Olap_vs_DMSO
require_median_log "K5_K19_Olap_vs_DMSO.log"
