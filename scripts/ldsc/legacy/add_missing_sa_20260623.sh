#!/bin/bash
# ============================================================
# Supplemental analysis: generate sumstats and calculate h2 and rg for two omitted targets.
# Write sumstats to a writable OUT_DIR and append results to the existing CSV.
# ============================================================
# Usage: nohup bash add_missing_targets.sh > add_missing.log 2>&1 &
# ============================================================

# ==================== Configuration ====================
TARGET_DIR="/exeh_3/rstao/dr.so/002.meta/002.clear.data/gc.targets"
LIBRARY_DIR="/exeh_3/rstao/dr.so/002.meta/002.clear.data/gc.library"
LDSC="/exeh_4/jinghong_qiu/ldsc/ldsc.py"
MUNGE="/exeh_4/jinghong_qiu/ldsc/munge_sumstats.py"
REF_LD="/exeh_3/rstao/programs/ldsc/eur_w_ld_chr/"
W_LD="/exeh_3/rstao/programs/ldsc/eur_w_ld_chr/"
OUT_DIR="/exeh_4/jinghong_qiu/SO_Lab/15.fdrreg.update/01.result/10.ldsc"
LOG_DIR="${OUT_DIR}/logs"

# Output directory for munge products (must be writable).
MUNGE_OUT_DIR="${OUT_DIR}/munge"

# Match the original munge step by aligning to the HapMap3 SNP list.
MERGE_ALLELES="/exeh_3/rstao/programs/ldsc/w_hm3.snplist"

H2_CSV="${OUT_DIR}/heritability_results.csv"
RG_CSV="${OUT_DIR}/genetic_correlation_results.csv"

mkdir -p "$LOG_DIR" "$MUNGE_OUT_DIR"

# ==================== Omitted source files ====================
# Add additional files here when needed.
MISSING_RAW=(
    "/exeh_3/rstao/dr.so/002.meta/002.clear.data/targets/sa.bpd.pgc.clear.txt"
    "/exeh_3/rstao/dr.so/002.meta/002.clear.data/targets/sa.scz.pgc.clear.txt"
)

# ==================== Step 0: munge into sumstats.gz ====================
echo "====== Step 0: munge omitted targets ======"

NEW_FILES=()
NEW_NAMES=()

for raw in "${MISSING_RAW[@]}"; do
    base="$(basename "$raw")"                       # sa.bpd.pgc.clear.txt
    name="${base%.clear.txt}"                       # sa.bpd.pgc
    out_prefix="${MUNGE_OUT_DIR}/${base}"           # Output prefix.
    sumstats="${out_prefix}.sumstats.gz"

    if [ -f "$sumstats" ]; then
        echo "  [SKIP munge] $sumstats already exists"
    else
        echo "  [munge] $base"

        # ---- Match the original logic: read a single N from the first row ----
        n_col=$(head -1 "$raw" | tr ' \t' '\n\n' | grep -nx "n" | head -1 | cut -d: -f1)
        if [ -z "$n_col" ]; then
            echo "  [ERROR] $base has no n column; skipping"
            continue
        fi
        nval=$(awk -v c="$n_col" 'NR==2{print $c}' "$raw")
        echo "        N = $nval (from the first row of the n column)"

        python "$MUNGE" \
            --sumstats "$raw" \
            --N-col n \
            --merge-alleles "$MERGE_ALLELES" \
            --out "$out_prefix"
    fi

    if [ ! -f "$sumstats" ]; then
        echo "  [ERROR] munge failed and did not generate $sumstats"
        continue
    fi

    NEW_FILES+=("$sumstats")
    NEW_NAMES+=("$name")
done

NNEW=${#NEW_FILES[@]}
if [ "$NNEW" -eq 0 ]; then
    echo "[FATAL] No new sumstats were generated; exiting"
    exit 1
fi

# ==================== Collect all existing files ====================
# New files are in MUNGE_OUT_DIR and do not conflict with the scan below.
EXIST_FILES=()
EXIST_NAMES=()

for f in "$TARGET_DIR"/*.clear.txt.sumstats.gz; do
    [ -f "$f" ] || continue
    EXIST_FILES+=("$f")
    EXIST_NAMES+=("$(basename "$f" .clear.txt.sumstats.gz)")
done

for f in "$LIBRARY_DIR"/*.sumstats.gz; do
    [ -f "$f" ] || continue
    EXIST_FILES+=("$f")
    EXIST_NAMES+=("$(basename "$f" .sumstats.gz)")
done

NEXIST=${#EXIST_FILES[@]}
echo ""
echo "============================================================"
echo "  new targets: $NNEW"
echo "  existing files: $NEXIST"
echo "  rg pairs to add: new x existing $((NNEW * NEXIST)) + new x new $((NNEW*(NNEW-1)/2))"
echo "============================================================"
echo ""

# ==================== Step 1: h2 for new files ====================
echo "====== Step 1: heritability for new targets ======"

for (( k=0; k<NNEW; k++ )); do
    file="${NEW_FILES[$k]}"
    name="${NEW_NAMES[$k]}"
    src="T"
    out="${LOG_DIR}/${src}__${name}"
    log="${out}.log"

    if [ -f "$log" ] && grep -q "Analysis finished" "$log" 2>/dev/null; then
        echo "  [SKIP] ${src}:${name}"
    else
        echo "  [h2]  ${src}:${name}"
        python "$LDSC" --h2 "$file" \
            --ref-ld-chr "$REF_LD" \
            --w-ld-chr "$W_LD" \
            --out "$out"
    fi

    if [ ! -f "$log" ]; then
        echo "  [ERROR] ${src}:${name} - log was not generated"
        echo "${name},${src},NA,NA,NA,NA" >> "$H2_CSV"
        continue
    fi

    lambda_gc=$(grep "Lambda GC:"    "$log" | awk '{print $NF}')
    mean_chi2=$(grep "Mean Chi^2:"   "$log" | awk '{print $NF}')
    intercept=$(grep "Intercept:"    "$log" | awk '{print $2}')
    intercept_se=$(grep "Intercept:" "$log" | sed 's/.*(\(.*\))/\1/')

    echo "${name},${src},${lambda_gc},${mean_chi2},${intercept},${intercept_se}" >> "$H2_CSV"
done

# ==================== Step 2: calculate missing rg ====================
echo ""
echo "====== Step 2: calculate missing Genetic Correlations ======"

parse_and_append_rg () {
    local log="$1" name_i="$2" name_j="$3"

    if [ ! -f "$log" ] || ! grep -q "Genetic Correlation:" "$log" 2>/dev/null; then
        echo "        [ERROR] failed; inspect: $log"
        echo "${name_i},${name_j},NA,NA,NA,NA,NA,NA" >> "$RG_CSV"
        return
    fi

    local rg se rg_pval gcov_line gcov_intercept gcov_intercept_se gcov_int_pval
    rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
    se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
    rg_pval=$(grep "^P:" "$log" | awk '{print $2}')

    gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
    gcov_intercept=$(echo "$gcov_line" | awk '{print $2}')
    gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/')

    gcov_int_pval=$(python -c "
from scipy.stats import norm
import sys
z = float(sys.argv[1]) / float(sys.argv[2])
print('%e' % (2.0 * norm.sf(abs(z))))
" "$gcov_intercept" "$gcov_intercept_se" 2>/dev/null || echo "NA")

    echo "${name_i},${name_j},${rg},${se},${rg_pval},${gcov_intercept},${gcov_intercept_se},${gcov_int_pval}" >> "$RG_CSV"
}

NPAIR=$(( NNEW * NEXIST + NNEW * (NNEW - 1) / 2 ))
count=0

# ---- (A) New files versus all existing files ----
for (( a=0; a<NNEW; a++ )); do
    name_i="${NEW_NAMES[$a]}"; src_i="T"
    for (( b=0; b<NEXIST; b++ )); do
        ((count++))
        name_j="${EXIST_NAMES[$b]}"
        case "${EXIST_FILES[$b]}" in
            "$TARGET_DIR"/*) src_j="T" ;;
            *)               src_j="L" ;;
        esac

        out="${LOG_DIR}/${src_i}__${name_i}__vs__${src_j}__${name_j}"
        log="${out}.log"
        echo "  [$count/$NPAIR] ${src_i}:${name_i} vs ${src_j}:${name_j}"

        if [ -f "$log" ] && grep -q "Analysis finished" "$log" 2>/dev/null; then
            echo "        [SKIP] already complete"
        else
            python "$LDSC" --rg "${NEW_FILES[$a]},${EXIST_FILES[$b]}" \
                --ref-ld-chr "$REF_LD" \
                --w-ld-chr "$W_LD" \
                --out "$out"
        fi
        parse_and_append_rg "$log" "$name_i" "$name_j"
    done
done

# ---- (B) All pairs among new files ----
for (( a=0; a<NNEW; a++ )); do
    for (( b=a+1; b<NNEW; b++ )); do
        ((count++))
        name_i="${NEW_NAMES[$a]}"; src_i="T"
        name_j="${NEW_NAMES[$b]}"; src_j="T"

        out="${LOG_DIR}/${src_i}__${name_i}__vs__${src_j}__${name_j}"
        log="${out}.log"
        echo "  [$count/$NPAIR] ${src_i}:${name_i} vs ${src_j}:${name_j}"

        if [ -f "$log" ] && grep -q "Analysis finished" "$log" 2>/dev/null; then
            echo "        [SKIP] already complete"
        else
            python "$LDSC" --rg "${NEW_FILES[$a]},${NEW_FILES[$b]}" \
                --ref-ld-chr "$REF_LD" \
                --w-ld-chr "$W_LD" \
                --out "$out"
        fi
        parse_and_append_rg "$log" "$name_i" "$name_j"
    done
done

echo ""
echo "============================================================"
echo "  Supplemental analysis complete"
echo "  sumstats output: $MUNGE_OUT_DIR/"
echo "  Heritability:  $H2_CSV"
echo "  Genetic Corr:  $RG_CSV"
echo "============================================================"
