#!/bin/bash
# ============================================================
# Test script: validate parsing and pairing logic.
# Run: bash test_ldsc.sh
# ============================================================

PASS=0; FAIL=0
check() {
    local desc="$1" expected="$2" actual="$3"
    expected=$(echo "$expected" | tr -d '[:space:]')
    actual=$(echo "$actual" | tr -d '[:space:]')
    if [ "$expected" = "$actual" ]; then
        echo "  ✓ $desc = $actual"
        ((PASS++))
    else
        echo "  ✗ $desc: expected='$expected', got='$actual'"
        ((FAIL++))
    fi
}

TEST_DIR=$(mktemp -d)

# ==================== Part 1: h2 log parsing ====================
echo "=== Part 1: h2 log parsing ==="

cat > "$TEST_DIR/phenA.h2.log" << 'EOF'
Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.4238 (0.0146)
Lambda GC: 1.6908
Mean Chi^2: 1.9417
Intercept: 1.0474 (0.0118)
Ratio: 0.0504 (0.0126)

Analysis finished at Tue Aug 25 16:47:25 2020
Total time elapsed: 15.39s
EOF

log="$TEST_DIR/phenA.h2.log"
check "Lambda GC"    "1.6908" "$(grep "Lambda GC:" "$log" | awk '{print $NF}')"
check "Mean Chi^2"   "1.9417" "$(grep "Mean Chi^2:" "$log" | awk '{print $NF}')"
check "h2 Intercept" "1.0474" "$(grep "Intercept:" "$log" | awk '{print $2}')"
check "h2 Int SE"    "0.0118" "$(grep "Intercept:" "$log" | sed 's/.*(\(.*\))/\1/')"

# ==================== Part 2: rg log parsing ====================
echo ""
echo "=== Part 2: rg log parsing ==="

cat > "$TEST_DIR/A__vs__B.rg.log" << 'EOF'
Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.4238 (0.0146)
Lambda GC: 1.6908
Mean Chi^2: 1.9417
Intercept: 1.0474 (0.0118)
Ratio: 0.0504 (0.0126)

Heritability of phenotype 2/2
-----------------------------
Total Observed scale h2: 0.0221 (0.0043)
Lambda GC: 1.0466
Mean Chi^2: 1.0624
Intercept: 1.0072 (0.0073)
Ratio: 0.1158 (0.1172)

Genetic Covariance
------------------
Total Observed scale gencov: 0.013 (0.0047)
Mean z1*z2: 0.0241
Intercept: -0.0019 (0.0064)

Genetic Correlation
-------------------
Genetic Correlation: 0.134 (0.0508)
Z-score: 2.638
P: 0.0083


Summary of Genetic Correlation Results
p1                                              p2     rg      se      z       p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
/test/a.sumstats.gz  /test/b.sumstats.gz  0.134  0.0508  2.638  0.0083  0.0221     0.0043  1.0072     0.0073   -0.0019       0.0064

Analysis finished at Tue Aug 25 16:47:25 2020
Total time elapsed: 15.39s
EOF

log="$TEST_DIR/A__vs__B.rg.log"
rg=$(grep "Genetic Correlation:" "$log" | head -1 | awk '{print $3}')
se=$(grep "Genetic Correlation:" "$log" | head -1 | sed 's/.*(\(.*\))/\1/')
rg_pval=$(grep "^P:" "$log" | awk '{print $2}')
gcov_line=$(awk '/^Genetic Covariance$/,/^Genetic Correlation$/{if(/Intercept:/) print}' "$log")
gcov_intercept=$(echo "$gcov_line" | awk '{print $2}')
gcov_intercept_se=$(echo "$gcov_line" | sed 's/.*(\(.*\))/\1/')

check "rg"              "0.134"   "$rg"
check "se"              "0.0508"  "$se"
check "rg p-value"      "0.0083"  "$rg_pval"
check "gcov intercept"  "-0.0019" "$gcov_intercept"
check "gcov int SE"     "0.0064"  "$gcov_intercept_se"

# intercept p-value
gcov_int_pval=$(python -c "
from scipy.stats import norm
z = -0.0019 / 0.0064
print(f'{2.0 * norm.sf(abs(z)):.6e}')
" 2>/dev/null)
echo "  [info] intercept p-value = $gcov_int_pval (requires scipy)"

# ==================== Part 3: pairing logic ====================
echo ""
echo "=== Part 3: pairing logic (including duplicate-name conflicts) ==="

mkdir -p "$TEST_DIR/gc.targets" "$TEST_DIR/gc.librarys"
touch "$TEST_DIR/gc.targets/bip.clear.txt.sumstats.gz"
touch "$TEST_DIR/gc.targets/scz.clear.txt.sumstats.gz"
touch "$TEST_DIR/gc.targets/mdd.clear.txt.sumstats.gz"
touch "$TEST_DIR/gc.librarys/bip.clear.txt.sumstats.gz"
touch "$TEST_DIR/gc.librarys/adhd.clear.txt.sumstats.gz"

TARGET_DIR="$TEST_DIR/gc.targets"
LIBRARY_DIR="$TEST_DIR/gc.librarys"

FILES=()
for f in "$TARGET_DIR"/*.clear.txt.sumstats.gz "$LIBRARY_DIR"/*.clear.txt.sumstats.gz; do
    [ -f "$f" ] && FILES+=("$f")
done

N=${#FILES[@]}
NP=$(( N * (N - 1) / 2 ))
check "total files" "5" "$N"
check "total pairs" "10" "$NP"

echo ""
echo "  All pairs:"
pair_count=0
for (( i=0; i<N; i++ )); do
    for (( j=i+1; j<N; j++ )); do
        ((pair_count++))
        ni=$(basename "${FILES[$i]}" .clear.txt.sumstats.gz)
        nj=$(basename "${FILES[$j]}" .clear.txt.sumstats.gz)
        # Label the source.
        [[ "${FILES[$i]}" == *"/gc.targets/"* ]] && ni="T:${ni}" || ni="L:${ni}"
        [[ "${FILES[$j]}" == *"/gc.targets/"* ]] && nj="T:${nj}" || nj="L:${nj}"
        echo "    $pair_count. $ni <-> $nj"
    done
done
check "actual pair count" "10" "$pair_count"

rm -rf "$TEST_DIR"

echo ""
echo "============================================"
echo "  Test results: PASS=$PASS, FAIL=$FAIL"
if [ "$FAIL" -eq 0 ]; then
    echo "  All tests passed"
else
    echo "  Some tests failed; review the output"
fi
echo "============================================"
