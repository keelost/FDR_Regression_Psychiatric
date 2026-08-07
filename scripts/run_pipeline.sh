#!/usr/bin/env bash
# Ordinary command-line runner for the prepared FDRreg analysis tree.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGE="all"
TARGET_SPEC="all"
VERSION="both"
JOBS=1
CONFIG_FILE="${ROOT_DIR}/config/project.env"

usage() {
  cat <<'EOF'
Usage: scripts/run_pipeline.sh [options]

Stages: imputation, ldsc, snp, magma, metaxcan, smultixcan, sensitivity, enrichment,
        validation, summary, all, list
Options:
  --stage NAME       Stage to run (default: all)
  --target LIST      Comma-separated targets, or all (default: all)
  --version VALUE    v7, v8, or both (default: both)
  --jobs N           Parallel jobs passed to versioned runners (default: 1)
  --config FILE      Environment configuration file
  -h, --help         Show this help

The imputation stage starts from cleaned <trait>.clear.txt files. Other stages
start from prepared per-target folders and generated *.impute.map.txt files.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage) [[ $# -ge 2 ]] || die "--stage needs a value"; STAGE="$2"; shift 2 ;;
    --target) [[ $# -ge 2 ]] || die "--target needs a value"; TARGET_SPEC="$2"; shift 2 ;;
    --version) [[ $# -ge 2 ]] || die "--version needs a value"; VERSION="$2"; shift 2 ;;
    --jobs) [[ $# -ge 2 ]] || die "--jobs needs a value"; JOBS="$2"; shift 2 ;;
    --config) [[ $# -ge 2 ]] || die "--config needs a value"; CONFIG_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -f "$CONFIG_FILE" ]] || [[ -f "${ROOT_DIR}/config/project.env.example" ]] || die "configuration file not found"
if [[ -f "$CONFIG_FILE" ]]; then source "$CONFIG_FILE"; else source "${ROOT_DIR}/config/project.env.example"; fi

FDRREG_BASE_DIR="${FDRREG_BASE_DIR:-${FDRREG_RESULTS_DIR:-}}"
FDRREG_RESULTS_DIR="${FDRREG_RESULTS_DIR:-$FDRREG_BASE_DIR}"
FDRREG_INPUT_DIR="${FDRREG_INPUT_DIR:-$FDRREG_BASE_DIR}"
FDRREG_R="${FDRREG_R:-Rscript}"
FDRREG_PYTHON="${FDRREG_PYTHON:-python}"
if [[ "$STAGE" != list ]]; then
  [[ -d "$FDRREG_BASE_DIR" ]] || die "FDRREG_BASE_DIR does not exist: $FDRREG_BASE_DIR"
fi
[[ "$VERSION" == v7 || "$VERSION" == v8 || "$VERSION" == both ]] || die "--version must be v7, v8, or both"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"

TARGET_FILE="${ROOT_DIR}/config/targets.tsv"
[[ -f "$TARGET_FILE" ]] || die "missing target configuration: $TARGET_FILE"

targets() {
  if [[ "$TARGET_SPEC" == all ]]; then awk -F '\t' 'NR > 1 {print $1}' "$TARGET_FILE"
  else tr ',' '\n' <<< "$TARGET_SPEC" | sed '/^[[:space:]]*$/d;s/^[[:space:]]*//;s/[[:space:]]*$//' ; fi
}

validate_targets() {
  local target
  while read -r target; do
    [[ -n "$target" ]] || continue
    target_field "$target" traits_with >/dev/null 2>&1 || die "unknown target '$target' (see config/targets.tsv)"
  done < <(targets)
}

target_field() {
  local target="$1" field="$2"
  awk -F '\t' -v t="$target" -v f="$field" 'NR == 1 { for (i=1;i<=NF;i++) if ($i==f) c=i; next } $1==t {print $c; found=1} END {if (!found) exit 1}' "$TARGET_FILE"
}

selected_versions() {
  case "$VERSION" in
    both) printf '%s\n' v8 v7 ;;
    v8|v7) printf '%s\n' "$VERSION" ;;
  esac
}

imputation_traits() {
  {
    targets
    while read -r target; do
      [[ -n "$target" ]] || continue
      target_field "$target" traits_no | tr ',' '\n'
      target_field "$target" traits_with | tr ',' '\n'
    done < <(targets)
  } | sed '/^[[:space:]]*$/d' | sort -u | paste -sd, -
}

run_imputation() {
  bash "$ROOT_DIR/scripts/imputation/run_dist_imputation.sh" \
    --traits "$(imputation_traits)" --jobs "$JOBS" --config "$CONFIG_FILE"
}

run_snp() {
  while read -r target; do
    [[ -n "$target" ]] || continue
    with="$(target_field "$target" traits_with || true)"
    without="$(target_field "$target" traits_no || true)"
    "$FDRREG_R" "$ROOT_DIR/scripts/real/common/fdrreg.R" \
      --target "$target" --traits_with "$with" --traits_no "$without" \
      --input_dir "$FDRREG_INPUT_DIR" --output_dir "$FDRREG_RESULTS_DIR" \
      --ldsc_results "${FDRREG_LDSC_RESULTS:-}" --seed "${FDRREG_SEED:-100}"
  done < <(targets)
}

run_magma() {
  [[ -n "${FDRREG_BIO_ENTREZ:-}" ]] || die "set FDRREG_BIO_ENTREZ for the MAGMA stage"
  while read -r target; do
    [[ -n "$target" ]] || continue
    with="$(target_field "$target" traits_with || true)"
    without="$(target_field "$target" traits_no || true)"
    vars="${with}${with:+${without:+,}}${without}"
    "$FDRREG_R" "$ROOT_DIR/scripts/real/common/magma_fdrreg_analysis.R" \
      --trait "$target" --variables "$vars" --bio "$FDRREG_BIO_ENTREZ" \
      --output "$FDRREG_RESULTS_DIR" --seed "${FDRREG_SEED:-100}"
  done < <(targets)
}

run_metaxcan() {
  [[ -n "${FDRREG_BIO_ENSEMBL:-}" ]] || die "set FDRREG_BIO_ENSEMBL for the MetaXcan stage"
  while read -r version; do
    [[ -n "$version" ]] || continue
    script="$ROOT_DIR/scripts/real/common/metaxcan_fdrreg.R"; [[ "$version" == v7 ]] && script="$ROOT_DIR/scripts/real/v7/metaxcan_fdrreg.R"
    "$FDRREG_R" "$script" --target "$TARGET_SPEC" --jobs "$JOBS" --input_dir "$FDRREG_RESULTS_DIR" \
      --output_dir "$FDRREG_RESULTS_DIR" --bio_file "$FDRREG_BIO_ENSEMBL"
  done < <(selected_versions)
}

run_smultixcan() {
  [[ -n "${FDRREG_BIO_ENSEMBL:-}" ]] || die "set FDRREG_BIO_ENSEMBL for the SMultiXcan stage"
  while read -r version; do
    [[ -n "$version" ]] || continue
    if [[ "$version" == v7 ]]; then script="$ROOT_DIR/scripts/real/v7/smultixcan_fdrreg_v7.R"; else script="$ROOT_DIR/scripts/real/common/smultixcan_fdrreg.R"; fi
    target_arg="$TARGET_SPEC"; [[ "$target_arg" == all ]] && target_arg="$(targets | paste -sd, -)"
    args=(--target "$target_arg" --base_dir "$FDRREG_RESULTS_DIR"
      --bio_library "$FDRREG_BIO_ENSEMBL" --seed "${FDRREG_SEED:-100}")
    [[ "$version" == v7 ]] && args+=(--jobs "$JOBS")
    "$FDRREG_R" "$script" "${args[@]}"
  done < <(selected_versions)
}

run_sensitivity() {
  while read -r target; do
    [[ -n "$target" ]] || continue
    with="$(target_field "$target" traits_with || true)"; without="$(target_field "$target" traits_no || true)"
    if [[ "$VERSION" == v8 || "$VERSION" == both ]]; then
      FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/common/fdrreg_sensitivity_analysis_rg.R" --target "$target" --traits_with "$with" --traits_no "$without"
      FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/common/fdrreg_sensitivity_analysis_z_abs.r" --target "$target" --traits_with "$with" --traits_no "$without"
    fi
    if [[ "$VERSION" == v7 || "$VERSION" == both ]]; then
      FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/v7/fdrreg_sensitivity_analysis_rg_v7.R" --target "$target" --traits_with "$with" --traits_no "$without"
      FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/v7/fdrreg_sensitivity_analysis_z_abs_v7.R" --target "$target" --traits_with "$with" --traits_no "$without"
    fi
  done < <(targets)
}

run_enrichment() {
  if [[ "$VERSION" == v8 || "$VERSION" == both ]]; then
    FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/common/DEA.R"
    FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/common/kegg_go_enrichment.r"
  fi
  if [[ "$VERSION" == v7 || "$VERSION" == both ]]; then
    FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/v7/DEA_v7.R"
    FDRREG_RESULTS_DIR="$FDRREG_RESULTS_DIR" "$FDRREG_R" "$ROOT_DIR/scripts/real/v7/kegg_go_enrichment_v7.r"
  fi
}

run_validation() {
  "$FDRREG_R" "$ROOT_DIR/scripts/real/common/validation_sensitivity_ppv.R" \
    --base-path "$FDRREG_RESULTS_DIR" \
    --base-output "$FDRREG_RESULTS_DIR/01.extra.analysis/02.sen_ppv"
}

run_summary() {
  "$FDRREG_PYTHON" "$ROOT_DIR/scripts/real/utilities/summarize_fdr_counts.py" \
    --root-dir "$FDRREG_RESULTS_DIR" \
    --risk-loci-dir "$FDRREG_RESULTS_DIR/01.extra.analysis/05.risk_loci" \
    --output-dir "$FDRREG_RESULTS_DIR/01.extra.analysis/09.summary_tables"
}

run_stage() {
  case "$1" in
    imputation) run_imputation ;;
    ldsc)
      bash "$ROOT_DIR/scripts/ldsc/run_ldsc_pipeline.sh"
      bash "$ROOT_DIR/scripts/ldsc/rebuild_results.sh"
      FDRREG_LDSC_RESULTS="${FDRREG_LDSC_OUT_DIR}/genetic_correlation_results.with_source.csv"
      ;;
    snp) run_snp ;;
    magma) run_magma ;;
    metaxcan) run_metaxcan ;;
    smultixcan) run_smultixcan ;;
    list) targets ;;
    sensitivity) run_sensitivity ;;
    enrichment) run_enrichment ;;
    validation) run_validation ;;
    summary) run_summary ;;
    all) run_snp; run_magma; run_metaxcan; run_smultixcan; run_sensitivity; run_enrichment; run_validation; run_summary ;;
    *) die "unknown stage: $1" ;;
  esac
}

validate_targets
run_stage "$STAGE"
