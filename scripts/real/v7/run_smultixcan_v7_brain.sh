#!/usr/bin/env bash
#
# Parallel batch S-MultiXcan (legacy MetaXcan v7) over GTEx v7 brain models.
#
# For each target, every trait that has S-PrediXcan v7 brain outputs under
# 06.metaxcan_v7/ is combined across the brain regions into one S-MultiXcan
# result. Jobs run in a parallel pool capped by -j.
#
# Usage:
#   run_smultixcan_v7_brain.sh                              # list targets and exit
#   run_smultixcan_v7_brain.sh --target mdd2019             # one target
#   run_smultixcan_v7_brain.sh -j 8 --target mdd2019,scz2014
#   run_smultixcan_v7_brain.sh -j 8 --target all
#
# The legacy SMulTiXcan.py runs under Python 2.7.15.
set -euo pipefail

# ----------------------------- Configuration -----------------------------
BASE_DIR="${FDRREG_RESULTS_DIR:?set FDRREG_RESULTS_DIR}"
SMULTIXCAN="${FDRREG_SMULTIXCAN:?set FDRREG_SMULTIXCAN}"
MODEL_DIR="${FDRREG_V7_MODEL_DIR:?set FDRREG_V7_MODEL_DIR}"
SNP_COV="${FDRREG_V7_SNP_COV:?set FDRREG_V7_SNP_COV}"

# v7 model file naming pattern (one capture group = tissue name).
MODELS_NAME_PATTERN='gtex_v7_(.*)_imputed_europeans_tw_0.5_signif.db'

# Python 2.7.15 interpreter for legacy MetaXcan. Used directly so parallel
# workers do not depend on `conda activate`.
# CHANGE THIS to your actual py2.7 env path.
PYTHON_BIN="${FDRREG_V7_PYTHON:-${FDRREG_PYTHON:-python}}"

# GWAS column names (must match the columns used in the S-PrediXcan v7 step).
SNP_COL="snpid"
EFFECT_ALLELE_COL="a1"
NON_EFFECT_ALLELE_COL="a2"
ZSCORE_COL="z.decor"
PVALUE_COL="p.decor"

# Expected number of brain-region S-PrediXcan files per trait.
N_BRAIN=13

# Default parallelism (override with -j).
JOBS=4
# --------------------------------------------------------------------------

usage() {
    echo "Usage: $(basename "$0") [-j N] --target t1[,t2,...] | --target all"
    echo "       $(basename "$0")            # list available targets and exit"
}

# ------------------------------ Argument parsing --------------------------
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
            TARGETS+=("$1"); shift ;;
    esac
done

# Discover every available target (directories that hold 06.metaxcan_v7).
ALL_TARGETS=()
for d in "${BASE_DIR}"/*/06.metaxcan_v7; do
    [[ -d "$d" ]] || continue
    ALL_TARGETS+=("$(basename "$(dirname "$d")")")
done

# No target -> list available targets and exit.
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    echo "Available targets under ${BASE_DIR}:"
    printf '  %s\n' "${ALL_TARGETS[@]}"
    echo
    echo "Run e.g.:  $(basename "$0") -j 8 --target mdd2019,scz2014"
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
[[ -f "$SMULTIXCAN" ]] || { echo "ERROR: SMulTiXcan not found: $SMULTIXCAN" >&2; exit 1; }
[[ -f "$SNP_COV" ]]    || { echo "ERROR: snp covariance not found: $SNP_COV" >&2; exit 1; }
[[ -x "$PYTHON_BIN" ]] || { echo "ERROR: python not found/executable: $PYTHON_BIN" >&2; exit 1; }
echo "Parallel jobs: ${JOBS}"

# --------------------------------------------------------------------------
# One S-MultiXcan job. One tab-separated arg:
#   metaxcan_folder <TAB> trait <TAB> gwas_file <TAB> out_file <TAB> log_file
smultixcan_one() {
    local IFS=$'\t'
    read -r mx_folder trait gwas out log <<< "$1"

    if [[ -s "$out" ]]; then
        echo "[${trait}] output exists, skip"
        return 0
    fi
    echo "[${trait}] running S-MultiXcan"
    "$PYTHON_BIN" "$SMULTIXCAN" \
        --models_folder "$MODEL_DIR/" \
        --models_name_pattern "$MODELS_NAME_PATTERN" \
        --snp_covariance "$SNP_COV" \
        --metaxcan_folder "$mx_folder" \
        --metaxcan_filter "gtex_v7_${trait}_in_Brain_.*\.csv" \
        --metaxcan_file_name_parse_pattern 'gtex_v7_(.*)_in_(.*).csv' \
        --gwas_file "$gwas" \
        --snp_column "$SNP_COL" \
        --effect_allele_column "$EFFECT_ALLELE_COL" \
        --non_effect_allele_column "$NON_EFFECT_ALLELE_COL" \
        --zscore_column "$ZSCORE_COL" \
        --pvalue_column "$PVALUE_COL" \
        --cutoff_condition_number 30 \
        --verbosity 7 \
        --throw \
        --output "$out" \
        >"$log" 2>&1 \
        && echo "[${trait}] done" \
        || { echo "[${trait}] FAILED (see $log)" >&2; return 1; }
}
export -f smultixcan_one
export PYTHON_BIN SMULTIXCAN MODEL_DIR SNP_COV MODELS_NAME_PATTERN \
       SNP_COL EFFECT_ALLELE_COL NON_EFFECT_ALLELE_COL ZSCORE_COL PVALUE_COL

# --------------------------------------------------------------------------
# Build job list across requested targets.
JOBS_LIST=()   # "mx_folder<TAB>trait<TAB>gwas<TAB>out<TAB>log"

for target in "${TARGETS[@]}"; do
    mx_folder="${BASE_DIR}/${target}/06.metaxcan_v7/"
    input_dir="${BASE_DIR}/${target}/01.magma_input"
    out_dir="${BASE_DIR}/${target}/08.smultixcan_output_v7"
    log_dir="${out_dir}/logs"
    if [[ ! -d "$mx_folder" ]]; then
        echo "[${target}] SKIP: no 06.metaxcan_v7" >&2
        continue
    fi
    mkdir -p "$out_dir" "$log_dir"

    # Discover traits from S-PrediXcan v7 filenames:
    #   gtex_v7_<trait>_in_Brain_<region>.csv
    shopt -s nullglob
    csvs=("${mx_folder}"gtex_v7_*_in_Brain_*.csv)
    shopt -u nullglob
    if [[ ${#csvs[@]} -eq 0 ]]; then
        echo "[${target}] SKIP: no brain S-PrediXcan v7 csv" >&2
        continue
    fi

    # Extract unique trait names.
    declare -A seen=()
    traits=()
    for f in "${csvs[@]}"; do
        b="$(basename "$f")"
        t="${b#gtex_v7_}"
        t="${t%%_in_Brain_*}"
        if [[ -z "${seen[$t]:-}" ]]; then
            seen[$t]=1
            traits+=("$t")
        fi
    done
    unset seen

    for trait in "${traits[@]}"; do
        # Count available brain csvs for this trait.
        shopt -s nullglob
        tfiles=("${mx_folder}gtex_v7_${trait}_in_Brain_"*.csv)
        shopt -u nullglob
        if [[ ${#tfiles[@]} -lt $N_BRAIN ]]; then
            echo "[${target}/${trait}] SKIP: only ${#tfiles[@]}/${N_BRAIN} brain csv present" >&2
            continue
        fi

        # v7 uses the original GWAS directly (no varID remap).
        gwas="${input_dir}/${trait}.overlap.4magma.txt"
        if [[ ! -s "$gwas" ]]; then
            echo "[${target}/${trait}] SKIP: GWAS missing: $gwas" >&2
            continue
        fi

        out="${out_dir}/${trait}.allbrain.txt"
        log="${log_dir}/${trait}.log"
        JOBS_LIST+=("${mx_folder}"$'\t'"${trait}"$'\t'"${gwas}"$'\t'"${out}"$'\t'"${log}")
    done
done

echo "==================== S-MultiXcan v7 (${#JOBS_LIST[@]} jobs) ===================="
if [[ ${#JOBS_LIST[@]} -gt 0 ]]; then
    printf '%s\0' "${JOBS_LIST[@]}" \
        | xargs -0 -P "$JOBS" -I{} bash -c 'smultixcan_one "$1"' _ {}
fi
echo "All requested targets finished."
