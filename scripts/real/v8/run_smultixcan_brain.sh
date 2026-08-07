#!/usr/bin/env bash
#
# Parallel batch S-MultiXcan over GTEx v8 mashr brain models (13 regions).
#
# For each target, every phenotype (trait) that has S-PrediXcan brain outputs
# under 06.metaxcan/ is combined across the 13 brain regions into one
# S-MultiXcan result. Jobs run in a parallel pool capped by -j.
#
# Usage:
#   run_smultixcan_brain.sh                       # list targets and exit
#   run_smultixcan_brain.sh mdd2019               # one target, default jobs
#   run_smultixcan_brain.sh -j 8 mdd2019 scz2014  # several targets, 8 parallel
#
set -euo pipefail

# ----------------------------- Configuration -----------------------------
BASE_DIR="${FDRREG_RESULTS_DIR:?set FDRREG_RESULTS_DIR}"
SMULTIXCAN="${FDRREG_SMULTIXCAN:?set FDRREG_SMULTIXCAN}"
MODEL_DIR="${FDRREG_V8_MODEL_DIR:?set FDRREG_V8_MODEL_DIR}"
SNP_COV="${FDRREG_V8_SNP_COV:?set FDRREG_V8_SNP_COV}"

# imlabtools interpreter (numpy 1.x). Used directly so parallel workers do
# not depend on `conda activate`.
PYTHON_BIN="${FDRREG_V8_PYTHON:-${FDRREG_PYTHON:-python}}"

# GWAS column names in the remapped *.varid.txt files
EFFECT_ALLELE_COL="a1"
NON_EFFECT_ALLELE_COL="a2"
ZSCORE_COL="z.decor"
PVALUE_COL="p.decor"

# Expected number of brain-region S-PrediXcan files per trait
N_BRAIN=13

# Default parallelism (override with -j)
JOBS=4
# --------------------------------------------------------------------------

# Parse -j.
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
    for d in "${BASE_DIR}"/*/06.metaxcan; do
        [[ -d "$d" ]] || continue
        echo "  $(basename "$(dirname "$d")")"
    done
    echo
    echo "Run e.g.:  $(basename "$0") -j 8 mdd2019 scz2014"
    exit 0
fi

# Sanity checks.
[[ -f "$SMULTIXCAN" ]] || { echo "ERROR: SMulTiXcan not found: $SMULTIXCAN" >&2; exit 1; }
[[ -f "$SNP_COV" ]]    || { echo "ERROR: snp covariance not found: $SNP_COV" >&2; exit 1; }
[[ -x "$PYTHON_BIN" ]] || { echo "ERROR: python not found: $PYTHON_BIN" >&2; exit 1; }
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
        --models_name_pattern 'mashr_(Brain_.*)\.db' \
        --model_db_snp_key varID \
        --snp_covariance "$SNP_COV" \
        --metaxcan_folder "$mx_folder" \
        --metaxcan_filter "gtex_v8_${trait}.overlap.4magma_in_Brain_.*\.csv" \
        --metaxcan_file_name_parse_pattern 'gtex_v8_(.*).overlap.4magma_in_(.*).csv' \
        --gwas_file "$gwas" \
        --separator $'\t' \
        --snp_column panel_variant_id \
        --effect_allele_column "$EFFECT_ALLELE_COL" \
        --non_effect_allele_column "$NON_EFFECT_ALLELE_COL" \
        --zscore_column "$ZSCORE_COL" \
        --pvalue_column "$PVALUE_COL" \
        --cutoff_condition_number 30 \
        --keep_non_rsid \
        --verbosity 7 \
        --throw \
        --output "$out" \
        >"$log" 2>&1 \
        && echo "[${trait}] done" \
        || { echo "[${trait}] FAILED (see $log)" >&2; return 1; }
}
export -f smultixcan_one
export PYTHON_BIN SMULTIXCAN MODEL_DIR SNP_COV \
       EFFECT_ALLELE_COL NON_EFFECT_ALLELE_COL ZSCORE_COL PVALUE_COL

# --------------------------------------------------------------------------
# Build job list across requested targets.
JOBS_LIST=()   # "mx_folder<TAB>trait<TAB>gwas<TAB>out<TAB>log"

for target in "$@"; do
    mx_folder="${BASE_DIR}/${target}/06.metaxcan/"
    out_dir="${BASE_DIR}/${target}/08.smultixcan_output"
    log_dir="${out_dir}/logs"
    if [[ ! -d "$mx_folder" ]]; then
        echo "[${target}] SKIP: no 06.metaxcan" >&2
        continue
    fi
    mkdir -p "$out_dir" "$log_dir"

    # Discover traits from S-PrediXcan filenames:
    #   gtex_v8_<trait>.overlap.4magma_in_Brain_<region>.csv
    shopt -s nullglob
    csvs=("${mx_folder}"gtex_v8_*.overlap.4magma_in_Brain_*.csv)
    shopt -u nullglob
    if [[ ${#csvs[@]} -eq 0 ]]; then
        echo "[${target}] SKIP: no brain S-PrediXcan csv" >&2
        continue
    fi

    # Extract unique trait names.
    declare -A seen=()
    traits=()
    for f in "${csvs[@]}"; do
        b="$(basename "$f")"
        t="${b#gtex_v8_}"
        t="${t%%.overlap.4magma_in_Brain_*}"
        if [[ -z "${seen[$t]:-}" ]]; then
            seen[$t]=1
            traits+=("$t")
        fi
    done
    unset seen

    for trait in "${traits[@]}"; do
        # Count available brain csvs for this trait.
        shopt -s nullglob
        tfiles=("${mx_folder}gtex_v8_${trait}.overlap.4magma_in_Brain_"*.csv)
        shopt -u nullglob
        if [[ ${#tfiles[@]} -lt $N_BRAIN ]]; then
            echo "[${target}/${trait}] SKIP: only ${#tfiles[@]}/${N_BRAIN} brain csv present" >&2
            continue
        fi

        gwas="${mx_folder}${trait}.overlap.4magma.varid.txt"
        if [[ ! -s "$gwas" ]]; then
            echo "[${target}/${trait}] SKIP: varid GWAS missing: $gwas" >&2
            continue
        fi

        out="${out_dir}/${trait}.allbrain.txt"
        log="${log_dir}/${trait}.log"
        JOBS_LIST+=("${mx_folder}"$'\t'"${trait}"$'\t'"${gwas}"$'\t'"${out}"$'\t'"${log}")
    done
done

echo "==================== S-MultiXcan (${#JOBS_LIST[@]} jobs) ===================="
if [[ ${#JOBS_LIST[@]} -gt 0 ]]; then
    printf '%s\0' "${JOBS_LIST[@]}" \
        | xargs -0 -P "$JOBS" -I{} bash -c 'smultixcan_one "$1"' _ {}
fi
echo "All requested targets finished."
