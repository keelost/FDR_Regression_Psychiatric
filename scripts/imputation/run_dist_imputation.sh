#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONFIG_FILE="${ROOT_DIR}/config/project.env"
TRAIT_SPEC=""
JOBS=1

usage() {
  cat <<'EOF'
Usage: scripts/imputation/run_dist_imputation.sh --traits LIST [options]

Options:
  --traits LIST      Comma-separated trait names (required)
  --jobs N           Maximum concurrent DIST chromosome jobs (default: 1)
  --config FILE      Environment configuration file
  -h, --help         Show this help

For each trait this creates a DIST input from <trait>.clear.txt, runs
chromosomes 1-22, maps non-rs identifiers, and writes
<FDRREG_INPUT_DIR>/<trait>.impute.map.txt.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --traits) [[ $# -ge 2 ]] || die "--traits needs a value"; TRAIT_SPEC="$2"; shift 2 ;;
    --jobs) [[ $# -ge 2 ]] || die "--jobs needs a value"; JOBS="$2"; shift 2 ;;
    --config) [[ $# -ge 2 ]] || die "--config needs a value"; CONFIG_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "$TRAIT_SPEC" ]] || die "--traits is required"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
[[ -f "$CONFIG_FILE" ]] || die "configuration file not found: $CONFIG_FILE"
source "$CONFIG_FILE"

R_BIN="${FDRREG_R:-Rscript}"
DIST_BIN="${FDRREG_DIST_BIN:-}"
CLEAR_TARGET_DIR="${FDRREG_CLEAR_TARGET_DIR:-}"
CLEAR_LIBRARY_DIR="${FDRREG_CLEAR_LIBRARY_DIR:-}"
WORK_DIR="${FDRREG_IMPUTE_WORK_DIR:-}"
OUTPUT_DIR="${FDRREG_INPUT_DIR:-}"
VARIANT_MAP="${FDRREG_DIST_VARIANT_MAP:-}"
REF_TABLE="${FDRREG_IMPUTATION_REFS:-${ROOT_DIR}/config/imputation_refs.tsv}"

[[ -x "$DIST_BIN" ]] || die "FDRREG_DIST_BIN is not executable: $DIST_BIN"
[[ -d "$CLEAR_TARGET_DIR" ]] || die "FDRREG_CLEAR_TARGET_DIR does not exist: $CLEAR_TARGET_DIR"
[[ -d "$CLEAR_LIBRARY_DIR" ]] || die "FDRREG_CLEAR_LIBRARY_DIR does not exist: $CLEAR_LIBRARY_DIR"
[[ -n "$WORK_DIR" ]] || die "set FDRREG_IMPUTE_WORK_DIR"
[[ -d "$OUTPUT_DIR" ]] || die "FDRREG_INPUT_DIR does not exist: $OUTPUT_DIR"
[[ -f "$VARIANT_MAP" ]] || die "FDRREG_DIST_VARIANT_MAP does not exist: $VARIANT_MAP"
[[ -f "$REF_TABLE" ]] || die "imputation reference table does not exist: $REF_TABLE"

reference_dir() {
  local trait="$1" panel var
  panel="$(awk -F '\t' -v t="$trait" 'NR>1 && $1==t {print $2; found=1} END{if(!found) exit 1}' "$REF_TABLE")" \
    || die "no imputation reference panel configured for trait '$trait'"
  var="FDRREG_DIST_REF_${panel}"
  printf '%s' "${!var:-}"
}

find_clear_file() {
  local trait="$1" candidate
  for candidate in \
    "$CLEAR_TARGET_DIR/${trait}.clear.txt" "$CLEAR_TARGET_DIR/${trait}.clear.txt.gz" \
    "$CLEAR_LIBRARY_DIR/${trait}.clear.txt" "$CLEAR_LIBRARY_DIR/${trait}.clear.txt.gz"; do
    [[ -f "$candidate" ]] && { printf '%s' "$candidate"; return 0; }
  done
  return 1
}

run_dist_chr() {
  local trait="$1" dist_input="$2" ref_dir="$3" chr="$4" panel="$5" pattern="$6" out ref_file
  out="$WORK_DIR/$trait/$trait.chr$chr.imputation.txt"
  ref_file="${pattern//\{chr\}/$chr}"
  ref_file="${ref_file//\{panel\}/${panel,,}}"
  [[ -f "$ref_dir/$ref_file" ]] || { echo "ERROR: DIST reference not found: $ref_dir/$ref_file" >&2; return 2; }
  "$DIST_BIN" "$dist_input" -c "$chr" -r "$ref_dir/$ref_file" -o "$out"
}
export DIST_BIN WORK_DIR
export -f run_dist_chr

IFS=',' read -r -a trait_list <<< "$TRAIT_SPEC"
for trait in "${trait_list[@]}"; do
  trait="${trait//[[:space:]]/}"
  [[ -n "$trait" ]] || continue
  clear_file="$(find_clear_file "$trait")" || die "cleaned GWAS not found for '$trait'"
  panel="$(awk -F '\t' -v t="$trait" 'NR>1 && $1==t {print $2}' "$REF_TABLE")"
  ref_dir="$(reference_dir "$trait")"
  pattern_var="FDRREG_DIST_REF_PATTERN_${panel}"
  ref_pattern="${!pattern_var:-chr{chr}.1kg.{panel}.gz}"
  [[ -d "$ref_dir" ]] || die "DIST $panel reference directory does not exist for '$trait': $ref_dir"
  mkdir -p "$WORK_DIR/$trait"
  dist_input="$WORK_DIR/$trait/$trait.4impute.txt"

  "$R_BIN" "$ROOT_DIR/scripts/imputation/prepare_dist_input.R" \
    --input "$clear_file" --output "$dist_input" \
    --clear-info-threshold "${FDRREG_CLEAR_INFO_THRESHOLD:-0.6}"

  printf '%s\n' {1..22} | xargs -P "$JOBS" -I CHR bash -c \
    'run_dist_chr "$1" "$2" "$3" "$4" "$5" "$6"' _ \
    "$trait" "$dist_input" "$ref_dir" CHR "$panel" "$ref_pattern"

  "$R_BIN" "$ROOT_DIR/scripts/imputation/merge_dist_output.R" \
    --trait "$trait" --input-dir "$WORK_DIR/$trait" --variant-map "$VARIANT_MAP" \
    --output "$OUTPUT_DIR/$trait.impute.map.txt" \
    --imputed-info-threshold "${FDRREG_IMPUTED_INFO_THRESHOLD:-0.5}"
done
