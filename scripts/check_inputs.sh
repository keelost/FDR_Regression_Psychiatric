#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${1:-${ROOT_DIR}/config/project.env}"
[[ -f "$CONFIG_FILE" ]] || CONFIG_FILE="${ROOT_DIR}/config/project.env.example"
source "$CONFIG_FILE"

fail=0
check_path() {
  local name="$1" value="${2:-}" kind="${3:-file}"
  if [[ -z "$value" ]]; then echo "MISSING  $name is not configured"; fail=1; return; fi
  if [[ "$kind" == dir && -d "$value" ]] || [[ "$kind" == file && -f "$value" ]]; then
    echo "OK       $name -> $value"
  else
    echo "MISSING  $name -> $value"
    fail=1
  fi
}

echo "Configured targets: $(awk -F '\t' 'NR>1{n++} END{print n+0}' "${ROOT_DIR}/config/targets.tsv")"
check_path FDRREG_BASE_DIR "${FDRREG_BASE_DIR:-}" dir
check_path FDRREG_INPUT_DIR "${FDRREG_INPUT_DIR:-}" dir
check_path FDRREG_LDSC_RESULTS "${FDRREG_LDSC_RESULTS:-}" file
check_path FDRREG_BIO_ENTREZ "${FDRREG_BIO_ENTREZ:-}" file
check_path FDRREG_BIO_ENSEMBL "${FDRREG_BIO_ENSEMBL:-}" file

if [[ -n "${FDRREG_DIST_BIN:-}${FDRREG_IMPUTE_WORK_DIR:-}" ]]; then
  check_path FDRREG_DIST_BIN "${FDRREG_DIST_BIN:-}" file
  check_path FDRREG_CLEAR_TARGET_DIR "${FDRREG_CLEAR_TARGET_DIR:-}" dir
  check_path FDRREG_CLEAR_LIBRARY_DIR "${FDRREG_CLEAR_LIBRARY_DIR:-}" dir
  check_path FDRREG_DIST_REF_EUR "${FDRREG_DIST_REF_EUR:-}" dir
  check_path FDRREG_DIST_REF_EAS "${FDRREG_DIST_REF_EAS:-}" dir
  check_path FDRREG_DIST_VARIANT_MAP "${FDRREG_DIST_VARIANT_MAP:-}" file
fi

exit "$fail"
