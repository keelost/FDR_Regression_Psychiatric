#!/usr/bin/env bash
#
# Parallel batch S-PrediXcan (legacy MetaXcan v7) over GTEx v7 brain models.
#
# Single phase: every independent (target, gwas, tissue) association goes
# into one parallel pool capped by -j. No varID remapping is needed for v7;
# the GWAS "snpid" column is used directly.
#
# Usage:
#   run_metaxcan_v7_brain.sh                                # list targets and exit
#   run_metaxcan_v7_brain.sh --target scz2014               # one target
#   run_metaxcan_v7_brain.sh -j 8 --target scz2014,bd2012   # several, 8 parallel
#   run_metaxcan_v7_brain.sh -j 8 --target all              # every target found
#
# The legacy MetaXcan.py runs under Python 2.7.15.
set -euo pipefail

# ----------------------------- Configuration -----------------------------
BASE_DIR="${FDRREG_RESULTS_DIR:?set FDRREG_RESULTS_DIR}"
METAXCAN="${FDRREG_METAXCAN:?set FDRREG_METAXCAN}"
MODEL_DIR="${FDRREG_V7_MODEL_DIR:?set FDRREG_V7_MODEL_DIR}"

# Python 2.7.15 interpreter for legacy MetaXcan. Used directly so parallel
# subshells do not depend on `conda activate`.
# CHANGE THIS to your actual py2.7 env path.
PYTHON_BIN="${FDRREG_V7_PYTHON:-${FDRREG_PYTHON:-python}}"

# GWAS column names (assumed identical across all *.4magma.txt files)
SNP_COL="snpid"
EFFECT_ALLELE_COL="a1"
NON_EFFECT_ALLELE_COL="a2"
ZSCORE_COL="z.decor"
PVALUE_COL="p.decor"

# Default parallelism (override with -j)
JOBS=4
# --------------------------------------------------------------------------

usage() {
    echo "Usage: $(basename "$0") [-j N] --target t1[,t2,...] | --target all"
    echo "       $(basename "$0")            # list available targets and exit"
}

# ------------------------------ Argument parsing --------------------------
# Manual parser so we can support the long option --target alongside -j.
TARGETS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            [[ $# -ge 2 ]] || { echo "Option $1 needs an argument" >&2; exit 1; }
            JOBS="$2"; shift 2 ;;
        --target)
            [[ $# -ge 2 ]] || { echo "Option --target needs an argument" >&2; exit 1; }
            IFS=',' read -r -a _t <<< "$2"
            TARGETS+=("${_t[@]}"); shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        -*)
            echo "Unknown option $1" >&2; usage; exit 1 ;;
        *)
            # Bare positional target, kept for backward compatibility.
            TARGETS+=("$1"); shift ;;
    esac
done

# Discover every available target (directories that hold 01.magma_input).
ALL_TARGETS=()
for d in "${BASE_DIR}"/*/01.magma_input; do
    [[ -d "$d" ]] || continue
    ALL_TARGETS+=("$(basename "$(dirname "$d")")")
done

# No target -> list available targets and exit.
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    echo "Available targets under ${BASE_DIR}:"
    printf '  %s\n' "${ALL_TARGETS[@]}"
    echo
    echo "Run e.g.:  $(basename "$0") -j 8 --target scz2014,bd2012"
    exit 0
fi

# Expand the special keyword "all".
for t in "${TARGETS[@]}"; do
    if [[ "$t" == "all" ]]; then
        TARGETS=("${ALL_TARGETS[@]}")
        break
    fi
done

# --------------------------------- Sanity checks --------------------------
[[ -f "$METAXCAN" ]]   || { echo "ERROR: MetaXcan.py not found: $METAXCAN" >&2; exit 1; }
[[ -x "$PYTHON_BIN" ]] || { echo "ERROR: python not found/executable: $PYTHON_BIN" >&2; exit 1; }

mapfile -t BRAIN_DBS < <(ls "${MODEL_DIR}"/gtex_v7_Brain_*_imputed_europeans_tw_0.5_signif.db 2>/dev/null | sort)
[[ ${#BRAIN_DBS[@]} -gt 0 ]] || { echo "ERROR: no gtex_v7_Brain_*.db in ${MODEL_DIR}" >&2; exit 1; }
echo "Found ${#BRAIN_DBS[@]} v7 brain models. Parallel jobs: ${JOBS}"

# --------------------------------------------------------------------------
# Worker for one association.
# One tab-separated arg: db <TAB> cov <TAB> gwas <TAB> out_csv <TAB> label
assoc_one() {
    local IFS=$'\t'
    read -r db cov gwas out_csv label <<< "$1"
    if [[ -s "$out_csv" ]]; then
        echo "[assoc] ${label}: exists, skip"
        return 0
    fi
    echo "[assoc] ${label}: running"
    "$PYTHON_BIN" "$METAXCAN" \
        --model_db_path "$db" \
        --covariance "$cov" \
        --gwas_file "$gwas" \
        --snp_column "$SNP_COL" \
        --effect_allele_column "$EFFECT_ALLELE_COL" \
        --non_effect_allele_column "$NON_EFFECT_ALLELE_COL" \
        --zscore_column "$ZSCORE_COL" \
        --pvalue_column "$PVALUE_COL" \
        --output_file "$out_csv" \
        >/dev/null 2>&1 \
        || { echo "[assoc] ${label}: FAILED" >&2; return 1; }
    echo "[assoc] ${label}: done"
}

export -f assoc_one
export PYTHON_BIN METAXCAN SNP_COL \
       EFFECT_ALLELE_COL NON_EFFECT_ALLELE_COL ZSCORE_COL PVALUE_COL

# --------------------------------------------------------------------------
# Build the full job list across all requested targets.
ASSOC_JOBS=()   # "db<TAB>cov<TAB>gwas<TAB>out_csv<TAB>label"

for target in "${TARGETS[@]}"; do
    input_dir="${BASE_DIR}/${target}/01.magma_input"
    out_dir="${BASE_DIR}/${target}/06.metaxcan_v7"
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
        # Short trait label: drop the ".overlap.4magma.txt" tail.
        gb="$(basename "$gwas")"
        gb="${gb%.txt}"
        gb="${gb%.overlap.4magma}"

        for db in "${BRAIN_DBS[@]}"; do
            # gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif -> Brain_Amygdala
            base="$(basename "$db" .db)"
            tissue="${base#gtex_v7_}"
            tissue="${tissue%_imputed_europeans_tw_0.5_signif}"
            cov="${MODEL_DIR}/gtex_v7_${tissue}_imputed_eur_covariances.txt.gz"
            out_csv="${out_dir}/gtex_v7_${gb}_in_${tissue}.csv"
            [[ -f "$cov" ]] || { echo "[${target}/${gb}/${tissue}] SKIP: cov missing" >&2; continue; }
            label="${target}/${gb}/${tissue}"
            ASSOC_JOBS+=("${db}"$'\t'"${cov}"$'\t'"${gwas}"$'\t'"${out_csv}"$'\t'"${label}")
        done
    done
done

# ------------------------------- Run associations -------------------------
echo "==================== Association (${#ASSOC_JOBS[@]} jobs) ===================="
if [[ ${#ASSOC_JOBS[@]} -gt 0 ]]; then
    printf '%s\0' "${ASSOC_JOBS[@]}" \
        | xargs -0 -P "$JOBS" -I{} bash -c 'assoc_one "$1"' _ {}
fi

echo "All requested targets finished."
