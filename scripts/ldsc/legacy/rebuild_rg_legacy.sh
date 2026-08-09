#!/bin/bash
# ============================================================
# Rebuild genetic_correlation_results.csv from existing logs.
# p1/p2 format examples: T:adhd2017 and L:bd2018.
# ============================================================

OUT_DIR="/path/to/SO_Lab/15.fdrreg.update/01.result/10.ldsc"
LOG_DIR="${OUT_DIR}/logs"
RG_CSV="${OUT_DIR}/genetic_correlation_results.csv"

echo "p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval" > "$RG_CSV"

declare -A seen
count=0; dup=0; skip=0

for log in "$LOG_DIR"/*__vs__*.log; do
    [ -f "$log" ] || continue
    fname=$(basename "$log" .log)

    # ---- Skip incomplete logs ----
    if ! grep -q "Analysis finished" "$log" 2>/dev/null; then
        ((skip++)); continue
    fi
    if ! grep -q "Genetic Correlation:" "$log" 2>/dev/null; then
        ((skip++)); continue
    fi

    # ---- Parse p1 and p2, including their sources ----
    left="${fname%%__vs__*}"
    right="${fname#*__vs__}"

    if [[ "$left" =~ ^[TL]__ ]] && [[ "$right" =~ ^[TL]__ ]]; then
        # New format: T__adhd2017__vs__L__bd2018.
        p1="${left:0:1}:${left:3}"
        p2="${right:0:1}:${right:3}"
    else
        # Old format: infer sources for adhd2017__vs__bd2018 from summary-table paths.
        summary_line=$(grep "\.sumstats" "$log" | tail -1)
        path1=$(echo "$summary_line" | awk '{print $1}')
        path2=$(echo "$summary_line" | awk '{print $2}')

        [[ "$path1" == *"/gc.targets/"* ]] && src1="T" || src1="L"
        [[ "$path2" == *"/gc.targets/"* ]] && src2="T" || src2="L"

        p1="${src1}:${left}"
        p2="${src2}:${right}"
    fi

    # ---- Deduplicate, preferring the newer log naming format ----
    key="${p1}|${p2}"
    if [ -n "${seen[$key]+x}" ]; then
        ((dup++)); continue
    fi
    seen[$key]=1
    ((count++))

    # ---- Parse rg ----
    rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
    se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
    rg_pval=$(grep "^P:" "$log" | awk '{print $2}')

    # ---- Parse the genetic-covariance intercept ----
    gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
    gcov_intercept=$(echo "$gcov_line" | awk '{print $2}')
    gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/')

    # ---- Intercept p-value (Python 2.7 compatible) ----
    gcov_int_pval=$(python -c "
from scipy.stats import norm
import sys
z = float(sys.argv[1]) / float(sys.argv[2])
print('%e' % (2.0 * norm.sf(abs(z))))
" "$gcov_intercept" "$gcov_intercept_se" 2>/dev/null || echo "NA")

    echo "${p1},${p2},${rg},${se},${rg_pval},${gcov_intercept},${gcov_intercept_se},${gcov_int_pval}" >> "$RG_CSV"
done

echo ""
echo "============================================================"
echo "  Complete"
echo "  pairs written: $count"
echo "  duplicates removed: $dup (old/new log overlap)"
echo "  skipped: $skip (incomplete or missing result)"
echo "  output: $RG_CSV"
echo "============================================================"
echo ""
echo "First 15 rows:"
head -15 "$RG_CSV"
