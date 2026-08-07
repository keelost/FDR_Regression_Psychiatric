#!/bin/bash

LDSC_OUT_DIR="${FDRREG_LDSC_OUT_DIR:?set FDRREG_LDSC_OUT_DIR}"
LOG_DIR="${LDSC_OUT_DIR}/logs"
OUT_CSV="${LDSC_OUT_DIR}/genetic_correlation_results.with_source.csv"
PYTHON_BIN="${FDRREG_LDSC_PYTHON:-${FDRREG_PYTHON:-python}}"

echo "p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval" > "$OUT_CSV"

shopt -s nullglob

count=0
for log in "$LOG_DIR"/*.log; do
  # 只处理 genetic correlation 结果对应的 log 文件，通常会有 "__vs__" 关键字
  if [[ "$log" != *"__vs__"* ]]; then
    continue
  fi

  filename=$(basename "$log" .log)
  # 命名形如 T__adhd2017__vs__L__bd2018
  # 解析出 p1_source, p1_name, p2_source, p2_name
  if [[ "$filename" =~ ^([TL])__([^_]+)__vs__([TL])__([^_]+)$ ]]; then
    p1_src="${BASH_REMATCH[1]}"
    p1_name="${BASH_REMATCH[2]}"
    p2_src="${BASH_REMATCH[3]}"
    p2_name="${BASH_REMATCH[4]}"
  else
    echo "[WARN] 跳过命名不符的文件 $log"
    continue
  fi

  # 检查日志是否含 Genetic Correlation 结果关键词
  if ! grep -q "Genetic Correlation:" "$log"; then
    echo "[WARN] $filename 无 Genetic Correlation 结果，跳过"
    continue
  fi

  # 解析 Genetic Correlation (寻找第一条)
  rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
  se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
  rg_pval=$(grep "^P:" "$log" | awk '{print $2}')

  gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
  gcov_intercept=$(echo "$gcov_line" | awk '{print $2}')
  gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/')

  gcov_int_pval=$("$PYTHON_BIN" -c "
from scipy.stats import norm
import sys
z = float(sys.argv[1]) / float(sys.argv[2])
print('%e' % (2.0 * norm.sf(abs(z))))
" "$gcov_intercept" "$gcov_intercept_se" 2>/dev/null || echo "NA")

  # 写入结果，格式化p1,p2为带源前缀格式
  echo "${p1_src}:${p1_name},${p2_src}:${p2_name},${rg},${se},${rg_pval},${gcov_intercept},${gcov_intercept_se},${gcov_int_pval}" >> "$OUT_CSV"

  ((count++))
done

echo "日志解析完成，共处理 $count 个 rg 结果，输出到 $OUT_CSV"
