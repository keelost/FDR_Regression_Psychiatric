#!/bin/bash
# ============================================================
# LDSC: heritability and all pairwise genetic correlations.
# ============================================================
# Usage: nohup bash run_ldsc_pipeline.sh > pipeline.log 2>&1 &
# ============================================================

# ==================== Configuration ====================
TARGET_DIR="${FDRREG_LDSC_TARGET_DIR:?set FDRREG_LDSC_TARGET_DIR}"
LIBRARY_DIR="${FDRREG_LDSC_LIBRARY_DIR:?set FDRREG_LDSC_LIBRARY_DIR}"
LDSC="${FDRREG_LDSC_BIN:?set FDRREG_LDSC_BIN to ldsc.py}"
REF_LD="${FDRREG_LDSC_REF:?set FDRREG_LDSC_REF}"
W_LD="${FDRREG_LDSC_WEIGHTS:-$REF_LD}"
OUT_DIR="${FDRREG_LDSC_OUT_DIR:?set FDRREG_LDSC_OUT_DIR}"
PYTHON_BIN="${FDRREG_LDSC_PYTHON:-${FDRREG_PYTHON:-python}}"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "$LOG_DIR"

# ==================== Collect all sumstats.gz files ====================
FILES=()
NAMES=()
SOURCES=()

# Targets use the suffix *.clear.txt.sumstats.gz.
for f in "$TARGET_DIR"/*.clear.txt.sumstats.gz; do
    [ -f "$f" ] || continue
    FILES+=("$f")
    NAMES+=("$(basename "$f" .clear.txt.sumstats.gz)")
    SOURCES+=("T")
done

# Library files use *.sumstats.gz (different from target names).
for f in "$LIBRARY_DIR"/*.sumstats.gz; do
    [ -f "$f" ] || continue
    FILES+=("$f")
    NAMES+=("$(basename "$f" .sumstats.gz)")
    SOURCES+=("L")
done

N=${#FILES[@]}
NT=$(printf '%s\n' "${SOURCES[@]}" | grep -c "^T$")
NL=$(printf '%s\n' "${SOURCES[@]}" | grep -c "^L$")
NP=$(( N * (N - 1) / 2 ))

echo "============================================================"
echo "  targets: $NT, library traits: $NL, total files: $N"
echo "  pairwise genetic correlations to calculate: $NP"
echo "============================================================"
echo ""

# ==================== CSV output ====================
H2_CSV="${OUT_DIR}/heritability_results.csv"
echo "file,source,lambda_GC,mean_chi2,h2_intercept,h2_intercept_se" > "$H2_CSV"

RG_CSV="${OUT_DIR}/genetic_correlation_results.csv"
echo "p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval" > "$RG_CSV"

# ==================== Step 1: Heritability ====================
echo "====== Step 1: Heritability ======"

for (( k=0; k<N; k++ )); do
    file="${FILES[$k]}"
    name="${NAMES[$k]}"
    src="${SOURCES[$k]}"
    unique="${src}__${name}"

    out="${LOG_DIR}/${unique}"
    log="${out}.log"

    if [ -f "$log" ] && grep -q "Analysis finished" "$log" 2>/dev/null; then
        echo "  [SKIP] ${src}:${name}"
    else
        echo "  [h2]  ${src}:${name}"
        "$PYTHON_BIN" "$LDSC" --h2 "$file" \
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

# ==================== Step 2: Pairwise Genetic Correlation ====================
echo ""
echo "====== Step 2: Pairwise Genetic Correlation ======"

count=0
for (( i=0; i<N; i++ )); do
    for (( j=i+1; j<N; j++ )); do
        ((count++))

        name_i="${NAMES[$i]}";  name_j="${NAMES[$j]}"
        src_i="${SOURCES[$i]}";  src_j="${SOURCES[$j]}"
        pair="${src_i}__${name_i}__vs__${src_j}__${name_j}"

        out="${LOG_DIR}/${pair}"
        log="${out}.log"

        echo "  [$count/$NP] ${src_i}:${name_i} vs ${src_j}:${name_j}"

        # ---- Skip completed jobs ----
        if [ -f "$log" ] && grep -q "Analysis finished" "$log" 2>/dev/null; then
            echo "        [SKIP] already complete"
        else
            "$PYTHON_BIN" "$LDSC" --rg "${FILES[$i]},${FILES[$j]}" \
                --ref-ld-chr "$REF_LD" \
                --w-ld-chr "$W_LD" \
                --out "$out"
        fi

        # ---- Check the log ----
        if [ ! -f "$log" ] || ! grep -q "Genetic Correlation:" "$log" 2>/dev/null; then
            echo "        [ERROR] failed; inspect: $log"
            echo "${name_i},${name_j},NA,NA,NA,NA,NA,NA" >> "$RG_CSV"
            continue
        fi

        # ---- Parse Genetic Correlation ----
        rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
        se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
        rg_pval=$(grep "^P:" "$log" | awk '{print $2}')

        # ---- Parse Genetic Covariance Intercept ----
        gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
        gcov_intercept=$(echo "$gcov_line" | awk '{print $2}')
        gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/')

        # ---- Intercept p-value (Python 2.7 compatible; arguments via sys.argv) ----
        gcov_int_pval=$(python -c "
from scipy.stats import norm
import sys
z = float(sys.argv[1]) / float(sys.argv[2])
print('%e' % (2.0 * norm.sf(abs(z))))
" "$gcov_intercept" "$gcov_intercept_se" 2>/dev/null || echo "NA")

        echo "${name_i},${name_j},${rg},${se},${rg_pval},${gcov_intercept},${gcov_intercept_se},${gcov_int_pval}" >> "$RG_CSV"
    done
done

echo ""
echo "============================================================"
echo "  All analyses complete"
echo "  Heritability:        $H2_CSV"
echo "  Genetic Correlation: $RG_CSV"
echo "  LDSC logs:           $LOG_DIR/"
echo "============================================================"
