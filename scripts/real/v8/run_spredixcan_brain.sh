#!/usr/bin/env bash
#
# Parallel batch S-PrediXcan over GTEx v8 mashr brain models (13 regions).
#
# Phase 1: remap each GWAS (rsID -> GTEx varID) once per file.
# Phase 2: run all (gwas x tissue) associations in a parallel pool.
#
# Parallelism spans both targets and tissues: every independent
# (target, gwas, tissue) job goes into one pool capped by -j.
#
# Usage:
#   run_spredixcan_brain.sh                       # list targets and exit
#   run_spredixcan_brain.sh scz2014               # one target, default jobs
#   run_spredixcan_brain.sh -j 8 scz2014 bd2012   # several targets, 8 parallel
# imlabtools conda env must be active for remap_gwas.py to work.
set -euo pipefail

# ----------------------------- Configuration -----------------------------
BASE_DIR="${FDRREG_RESULTS_DIR:?set FDRREG_RESULTS_DIR}"
METAXCAN="${FDRREG_METAXCAN:?set FDRREG_METAXCAN}"
MODEL_DIR="${FDRREG_V8_MODEL_DIR:?set FDRREG_V8_MODEL_DIR}"
CODE_DIR="${FDRREG_CODE_DIR:-$(cd "$(dirname "$0")/../utilities" && pwd)}"
REMAP_SCRIPT="${FDRREG_REMAP_SCRIPT:-${CODE_DIR}/remap_gwas.py}"
MAP_FILE="${FDRREG_MAP_FILE:-${CODE_DIR}/rsid_to_varid.tsv}"

# Python interpreter from the imlabtools env (numpy 1.x). Used directly so
# parallel subshells do not depend on `conda activate`.
PYTHON_BIN="${FDRREG_V8_PYTHON:-${FDRREG_PYTHON:-python}}"

# GWAS column names (assumed identical across all *.4magma.txt files)
SNP_COL="snpid"
EFFECT_ALLELE_COL="a1"
NON_EFFECT_ALLELE_COL="a2"
ZSCORE_COL="z.decor"
PVALUE_COL="p.decor"

# Default parallelism (override with -j)
JOBS=4
# --------------------------------------------------------------------------

# Parse -j flag.
while getopts ":j:" opt; do
    case "$opt" in
        j) JOBS="$OPTARG" ;;
        \?) echo "Unknown option -$OPTARG" >&2; exit 1 ;;
        :)  echo "Option -$OPTARG needs an argument" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# No target -> list available targets and exit.
if [[ $# -eq 0 ]]; then
    echo "Available targets under ${BASE_DIR}:"
    for d in "${BASE_DIR}"/*/01.magma_input; do
        [[ -d "$d" ]] || continue
        echo "  $(basename "$(dirname "$d")")"
    done
    echo
    echo "Run e.g.:  $(basename "$0") -j 8 scz2014 bd2012"
    exit 0
fi

# Sanity checks.
[[ -f "$MAP_FILE" ]]     || { echo "ERROR: map file not found: $MAP_FILE" >&2; exit 1; }
[[ -f "$REMAP_SCRIPT" ]] || { echo "ERROR: remap script not found: $REMAP_SCRIPT" >&2; exit 1; }
[[ -f "$METAXCAN" ]]     || { echo "ERROR: SPrediXcan not found: $METAXCAN" >&2; exit 1; }
[[ -x "$PYTHON_BIN" ]]   || { echo "ERROR: python not found: $PYTHON_BIN" >&2; exit 1; }

mapfile -t BRAIN_DBS < <(ls "${MODEL_DIR}"/mashr_Brain_*.db | sort)
[[ ${#BRAIN_DBS[@]} -gt 0 ]] || { echo "ERROR: no mashr_Brain_*.db in ${MODEL_DIR}" >&2; exit 1; }
echo "Found ${#BRAIN_DBS[@]} brain models. Parallel jobs: ${JOBS}"

# --------------------------------------------------------------------------
# Worker for one remap. Args: gwas_file  varid_out
remap_one() {
    local gwas="$1" out="$2"
    if [[ -s "$out" ]]; then
        echo "[remap] reuse $(basename "$out")"
        return 0
    fi
    echo "[remap] $(basename "$gwas") -> $(basename "$out")"
    "$PYTHON_BIN" "$REMAP_SCRIPT" \
        --gwas "$gwas" --map "$MAP_FILE" --snp_column "$SNP_COL" --out "$out" \
        >/dev/null 2>&1 \
        || { echo "[remap] FAILED: $gwas" >&2; return 1; }
}

# Worker for one association. One tab-separated arg: db <TAB> cov <TAB> varid_gwas <TAB> out_csv <TAB> label
assoc_one() {
    local IFS=$'\t'
    read -r db cov varid_gwas out_csv label <<< "$1"
    if [[ -s "$out_csv" ]]; then
        echo "[assoc] ${label}: exists, skip"
        return 0
    fi
    echo "[assoc] ${label}: running"
    "$PYTHON_BIN" "$METAXCAN" \
        --model_db_path "$db" \
        --model_db_snp_key varID \
        --covariance "$cov" \
        --gwas_file "$varid_gwas" \
        --separator $'\t' \
        --snp_column panel_variant_id \
        --effect_allele_column "$EFFECT_ALLELE_COL" \
        --non_effect_allele_column "$NON_EFFECT_ALLELE_COL" \
        --zscore_column "$ZSCORE_COL" \
        --pvalue_column "$PVALUE_COL" \
        --keep_non_rsid \
        --output_file "$out_csv" \
        >/dev/null 2>&1 \
        || { echo "[assoc] ${label}: FAILED" >&2; return 1; }
    echo "[assoc] ${label}: done"
}

export -f remap_one assoc_one
export PYTHON_BIN REMAP_SCRIPT MAP_FILE METAXCAN SNP_COL \
       EFFECT_ALLELE_COL NON_EFFECT_ALLELE_COL ZSCORE_COL PVALUE_COL

# --------------------------------------------------------------------------
# Build the full job list across all requested targets.
REMAP_JOBS=()   # "gwas<TAB>varid_out"
ASSOC_JOBS=()   # "db<TAB>cov<TAB>varid_gwas<TAB>out_csv<TAB>label"

for target in "$@"; do
    input_dir="${BASE_DIR}/${target}/01.magma_input"
    out_dir="${BASE_DIR}/${target}/06.metaxcan"
    if [[ ! -d "$input_dir" ]]; then
        echo "[${target}] SKIP: no 01.magma_input" >&2
        continue
    fi
    mkdir -p "$out_dir"

    shopt -s nullglob
    gwas_files=("${input_dir}"/*)
    shopt -u nullglob
    if [[ ${#gwas_files[@]} -eq 0 ]]; then
        echo "[${target}] SKIP: no GWAS files" >&2
        continue
    fi

    for gwas in "${gwas_files[@]}"; do
        gb="$(basename "$gwas")"; gb="${gb%.txt}"
        varid_gwas="${out_dir}/${gb}.varid.txt"
        REMAP_JOBS+=("${gwas}"$'\t'"${varid_gwas}")

        for db in "${BRAIN_DBS[@]}"; do
            tissue="$(basename "$db" .db)"; tissue="${tissue#mashr_}"
            cov="${MODEL_DIR}/$(basename "$db" .db).txt.gz"
            out_csv="${out_dir}/gtex_v8_${gb}_in_${tissue}.csv"
            [[ -f "$cov" ]] || { echo "[${target}/${gb}/${tissue}] SKIP: cov missing" >&2; continue; }
            label="${target}/${gb}/${tissue}"
            ASSOC_JOBS+=("${db}"$'\t'"${cov}"$'\t'"${varid_gwas}"$'\t'"${out_csv}"$'\t'"${label}")
        done
    done
done

# --------------------------- Phase 1: remap ------------------------------
echo "==================== Phase 1: remap (${#REMAP_JOBS[@]} files) ===================="
if [[ ${#REMAP_JOBS[@]} -gt 0 ]]; then
    printf '%s\0' "${REMAP_JOBS[@]}" \
        | xargs -0 -P "$JOBS" -I{} bash -c '
            IFS=$'"'"'\t'"'"' read -r g o <<< "$1"; remap_one "$g" "$o"' _ {}
fi

# ------------------------- Phase 2: association --------------------------
echo "==================== Phase 2: association (${#ASSOC_JOBS[@]} jobs) ===================="
if [[ ${#ASSOC_JOBS[@]} -gt 0 ]]; then
    printf '%s\0' "${ASSOC_JOBS[@]}" \
        | xargs -0 -P "$JOBS" -I{} bash -c 'assoc_one "$1"' _ {}
fi

echo "All requested targets finished."
