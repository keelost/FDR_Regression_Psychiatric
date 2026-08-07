#!/bin/bash
# ============================================================
# 重新解析已有 LDSC 日志 → 生成 CSV (无需重跑 LDSC)
# 运行: bash reparse_ldsc_logs.sh
# ============================================================

OUT_DIR="/exeh_4/jinghong_qiu/SO_Lab/15.fdrreg.update/01.result/10.ldsc"
LOG_DIR="${OUT_DIR}/logs"

# ==================== CSV 输出 ====================
H2_CSV="${OUT_DIR}/heritability_results.csv"
echo "file,source,lambda_GC,mean_chi2,h2_intercept,h2_intercept_se" > "$H2_CSV"

RG_CSV="${OUT_DIR}/genetic_correlation_results.csv"
echo "p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval" > "$RG_CSV"

h2_count=0
rg_count=0

# ==================== 遍历所有日志 ====================
for log in "$LOG_DIR"/*.log; do
    [ -f "$log" ] || continue
    log_name=$(basename "$log" .log)

    # ---------- 判断类型 ----------
    if grep -q "Genetic Correlation:" "$log" 2>/dev/null; then
        # ====== rg 日志 ======
        ((rg_count++))
        name_i="${log_name%%__vs__*}"
        name_j="${log_name##*__vs__}"

        if ! grep -q "Analysis finished" "$log" 2>/dev/null; then
            echo "  [rg] $name_i vs $name_j → 未完成, 跳过"
            echo "${name_i},${name_j},NA,NA,NA,NA,NA,NA" >> "$RG_CSV"
            continue
        fi

        # 解析 Genetic Correlation
        rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
        se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
        rg_pval=$(grep "^P:" "$log" | awk '{print $2}')

        # 解析 Genetic Covariance Intercept
        gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
        gcov_intercept=$(echo "$gcov_line" | awk '{print $2}' | tr -d '[:space:]')
        gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/' | tr -d '[:space:]')

        # 计算 intercept p-value (双侧检验) —— 修复版
        if [ -n "$gcov_intercept" ] && [ -n "$gcov_intercept_se" ]; then
            gcov_int_pval=$(python -c "from scipy.stats import norm; z=$gcov_intercept/$gcov_intercept_se; print('{:.6e}'.format(2.0*norm.sf(abs(z))))" 2>/dev/null)
            [ -z "$gcov_int_pval" ] && gcov_int_pval="NA"
        else
            gcov_int_pval="NA"
        fi

        echo "  [rg] $name_i vs $name_j: rg=$rg, p=$rg_pval, int=$gcov_intercept, int_p=$gcov_int_pval"
        echo "${name_i},${name_j},${rg},${se},${rg_pval},${gcov_intercept},${gcov_intercept_se},${gcov_int_pval}" >> "$RG_CSV"

    elif grep -q "Heritability of phenotype" "$log" 2>/dev/null; then
        # ====== h2 日志 ======
        ((h2_count++))

        # 判断来源和文件名
        if [[ "$log_name" == T__* ]]; then
            src="target";  fname="${log_name#T__}"
        elif [[ "$log_name" == L__* ]]; then
            src="library"; fname="${log_name#L__}"
        else
            src="unknown"; fname="$log_name"
        fi

        if ! grep -q "Analysis finished" "$log" 2>/dev/null; then
            echo "  [h2] $src/$fname → 未完成, 跳过"
            echo "${fname},${src},NA,NA,NA,NA" >> "$H2_CSV"
            continue
        fi

        lambda_gc=$(grep "Lambda GC:"   "$log" | awk '{print $NF}')
        mean_chi2=$(grep "Mean Chi^2:"  "$log" | awk '{print $NF}')
        intercept=$(grep "Intercept:"   "$log" | awk '{print $2}')
        intercept_se=$(grep "Intercept:" "$log" | sed 's/.*(\(.*\))/\1/')

        echo "  [h2] $src/$fname: lambda=$lambda_gc, intercept=$intercept ($intercept_se)"
        echo "${fname},${src},${lambda_gc},${mean_chi2},${intercept},${intercept_se}" >> "$H2_CSV"
    fi
done

echo ""
echo "============================================================"
echo "  重新解析完成!"
echo "  解析 h2 日志: $h2_count 个"
echo "  解析 rg 日志: $rg_count 个"
echo ""
echo "  Heritability:        $H2_CSV"
echo "  Genetic Correlation: $RG_CSV"
echo "============================================================"
