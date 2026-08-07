#!/bin/bash

LDSC_OUT_DIR="${FDRREG_LDSC_OUT_DIR:?set FDRREG_LDSC_OUT_DIR}"
LOG_DIR="${LDSC_OUT_DIR}/logs"
OUT_CSV="${LDSC_OUT_DIR}/genetic_correlation_results.with_source.csv"
PYTHON_BIN="${FDRREG_LDSC_PYTHON:-${FDRREG_PYTHON:-python}}"

echo "p1,p2,rg,se,rg_pval,gcov_intercept,gcov_intercept_se,gcov_intercept_pval" > "$OUT_CSV"

shopt -s nullglob

count=0
for log in "$LOG_DIR"/*.log; do
  # Process only genetic-correlation logs, normally identified by "__vs__".
  if [[ "$log" != *"__vs__"* ]]; then
    continue
  fi

  filename=$(basename "$log" .log)
  # Expected name: T__adhd2017__vs__L__bd2018.
  # Parse p1_source, p1_name, p2_source, and p2_name.
  if [[ "$filename" =~ ^([TL])__([^_]+)__vs__([TL])__([^_]+)$ ]]; then
    p1_src="${BASH_REMATCH[1]}"
    p1_name="${BASH_REMATCH[2]}"
    p2_src="${BASH_REMATCH[3]}"
    p2_name="${BASH_REMATCH[4]}"
  else
    echo "[WARN] Skipping file with an unexpected name: $log"
    continue
  fi

  # Check whether the log contains a Genetic Correlation result.
  if ! grep -q "Genetic Correlation:" "$log"; then
    echo "[WARN] $filename has no Genetic Correlation result; skipping"
    continue
  fi

  # Parse the first Genetic Correlation result.
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

  # Write the result with source-prefixed p1 and p2 identifiers.
  echo "${p1_src}:${p1_name},${p2_src}:${p2_name},${rg},${se},${rg_pval},${gcov_intercept},${gcov_intercept_se},${gcov_int_pval}" >> "$OUT_CSV"

  ((count++))
done

echo "Log parsing complete: processed $count rg results; output: $OUT_CSV"
