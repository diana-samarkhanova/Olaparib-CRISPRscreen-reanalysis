# Sourced by the RRA scripts; run from a complete checkout of scripts/mageck/.
# MAGeCK 0.5.9.5 can replace requested median normalization with total.
require_median_log() {
  local log_file="$1"
  if [[ ! -s "$log_file" ]]; then
    printf 'Cannot verify median normalization: missing log %s\n' "$log_file" >&2
    return 1
  fi
  if grep -Fq 'Switch to total read count normalization.' "$log_file"; then
    printf 'Rejected: MAGeCK switched to total normalization. Review %s; do not use this run as a median-normalized result.\n' "$log_file" >&2
    return 1
  fi
  if ! grep -Fq 'Final size factor:' "$log_file"; then
    printf 'Cannot verify normalization: no final size factors in %s\n' "$log_file" >&2
    return 1
  fi
}
