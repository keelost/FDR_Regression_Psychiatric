#!/bin/bash
# ============================================================
# 从已有日志重建 genetic_correlation_results.csv
# p1/p2 格式: T:adhd2017, L:bd2018 等
# ============================================================

OUT_DIR="/exeh_4/jinghong_qiu/SO_Lab/15.fdrreg.update/01.result/10.ldsc"
LOG_DIR="${OUT_DIR}/logs"
RG_CSV="${OUT_DIR}/genetic_correlation_results.csv"

echo "p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval" > "$RG_CSV"

declare -A seen
count=0; dup=0; skip=0

for log in "$LOG_DIR"/*__vs__*.log; do
    [ -f "$log" ] || continue
    fname=$(basename "$log" .log)

    # ---- 跳过不完整的日志 ----
    if ! grep -q "Analysis finished" "$log" 2>/dev/null; then
        ((skip++)); continue
    fi
    if ! grep -q "Genetic Correlation:" "$log" 2>/dev/null; then
        ((skip++)); continue
    fi

    # ---- 解析 p1, p2 (含来源) ----
    left="${fname%%__vs__*}"
    right="${fname#*__vs__}"

    if [[ "$left" =~ ^[TL]__ ]] && [[ "$right" =~ ^[TL]__ ]]; then
        # 新格式: T__adhd2017__vs__L__bd2018
        p1="${left:0:1}:${left:3}"
        p2="${right:0:1}:${right:3}"
    else
        # 旧格式: adhd2017__vs__bd2018 → 从日志 summary table 的路径推断
        summary_line=$(grep "\.sumstats" "$log" | tail -1)
        path1=$(echo "$summary_line" | awk '{print $1}')
        path2=$(echo "$summary_line" | awk '{print $2}')

        [[ "$path1" == *"/gc.targets/"* ]] && src1="T" || src1="L"
        [[ "$path2" == *"/gc.targets/"* ]] && src2="T" || src2="L"

        p1="${src1}:${left}"
        p2="${src2}:${right}"
    fi

    # ---- 去重 (新日志文件名靠前, 优先保留) ----
    key="${p1}|${p2}"
    if [ -n "${seen[$key]+x}" ]; then
        ((dup++)); continue
    fi
    seen[$key]=1
    ((count++))

    # ---- 解析 rg ----
    rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
    se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
    rg_pval=$(grep "^P:" "$log" | awk '{print $2}')

    # ---- 解析 gcov intercept ----
    gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
    gcov_intercept=$(echo "$gcov_line" | awk '{print $2}')
    gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/')

    # ---- intercept p-value (Python 2.7 兼容) ----
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
echo "  完成!"
echo "  写入: $count 对"
echo "  去重: $dup 个 (新旧日志重复)"
echo "  跳过: $skip 个 (未完成或无结果)"
echo "  输出: $RG_CSV"
echo "============================================================"
echo ""
echo "预览前 15 行:"
head -15 "$RG_CSV"
