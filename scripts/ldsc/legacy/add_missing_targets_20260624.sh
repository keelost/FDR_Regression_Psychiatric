#!/bin/bash
# ============================================================
# 补充分析: 为漏掉的 target 生成 sumstats 并补算 h2 + rg
# 支持两种 N 来源:
#   1) auto       —— 从原文件的 n 列读取 (与最初 munge 一致)
#   2) <整数>     —— 文件没有 n 列时, 使用固定样本量
# sumstats 输出到有写权限的 OUT_DIR, 结果追加到已有 CSV
# ============================================================
# Usage: nohup bash add_missing_targets.sh > add_missing.log 2>&1 &
# ============================================================

# ==================== 配置 ====================
TARGET_DIR="/exeh_3/rstao/dr.so/002.meta/002.clear.data/gc.targets"
LIBRARY_DIR="/exeh_3/rstao/dr.so/002.meta/002.clear.data/gc.library"
LDSC="/exeh_4/jinghong_qiu/ldsc/ldsc.py"
MUNGE="/exeh_4/jinghong_qiu/ldsc/munge_sumstats.py"
REF_LD="/exeh_3/rstao/programs/ldsc/eur_w_ld_chr/"
W_LD="/exeh_3/rstao/programs/ldsc/eur_w_ld_chr/"
OUT_DIR="/exeh_4/jinghong_qiu/SO_Lab/15.fdrreg.update/01.result/10.ldsc"
LOG_DIR="${OUT_DIR}/logs"

# munge 产物输出目录 (你有写权限的位置)
MUNGE_OUT_DIR="${OUT_DIR}/munge"

# 与你最初 munge 一致: 对齐 HapMap3 SNP 列表
MERGE_ALLELES="/exeh_3/rstao/programs/ldsc/w_hm3.snplist"

H2_CSV="${OUT_DIR}/heritability_results.csv"
RG_CSV="${OUT_DIR}/genetic_correlation_results.csv"

mkdir -p "$LOG_DIR" "$MUNGE_OUT_DIR"

# ==================== 漏掉的原始文件 ====================
# 格式: "<原始文件路径>|<N 来源>"
#   N 来源 = "auto"  —— 从文件第一行的 n 列自动取值
#   N 来源 = <整数>  —— 文件没有 n 列时, 直接指定样本量
# 已经 munge 过的文件会自动 SKIP, 重复列出无副作用
MISSING_RAW=(
    "/exeh_3/rstao/dr.so/002.meta/002.clear.data/targets/sa.bpd.pgc.clear.txt|auto"
    "/exeh_3/rstao/dr.so/002.meta/002.clear.data/targets/sa.scz.pgc.clear.txt|auto"
    "/exeh_3/rstao/dr.so/002.meta/002.clear.data/targets/adhd2016.clear.txt|17666"
)

# ==================== Step 0: munge 成 sumstats.gz ====================
echo "====== Step 0: munge 漏掉的 target ======"

NEW_FILES=()
NEW_NAMES=()

for entry in "${MISSING_RAW[@]}"; do
    # 解析 "path|N_spec"
    raw="${entry%|*}"
    nspec="${entry##*|}"

    base="$(basename "$raw")"                       # e.g. adhd2016.clear.txt
    name="${base%.clear.txt}"                       # e.g. adhd2016
    out_prefix="${MUNGE_OUT_DIR}/${base}"           # 输出前缀
    sumstats="${out_prefix}.sumstats.gz"

    if [ -f "$sumstats" ]; then
        echo "  [SKIP munge] $sumstats 已存在"
    else
        if [ ! -f "$raw" ]; then
            echo "  [ERROR] 原始文件不存在: $raw, 跳过"
            continue
        fi

        echo "  [munge] $base  (N spec: $nspec)"

        if [ "$nspec" = "auto" ]; then
            # ---- 与最初 munge 一致: 从 n 列第一行取单个 N ----
            n_col=$(head -1 "$raw" | tr ' \t' '\n\n' | grep -nx "n" | head -1 | cut -d: -f1)
            if [ -z "$n_col" ]; then
                echo "  [ERROR] $base 找不到 n 列, 跳过 (如该文件确无 n 列, 请把 N spec 改成具体整数)"
                continue
            fi
            nval=$(awk -v c="$n_col" 'NR==2{print $c}' "$raw")
            echo "        N = $nval (取自 n 列第一行)"

            python "$MUNGE" \
                --sumstats "$raw" \
                --N-col n \
                --merge-alleles "$MERGE_ALLELES" \
                --out "$out_prefix"
        else
            # ---- 文件无 n 列, 使用固定 N ----
            # 简单校验: 必须是正整数
            if ! [[ "$nspec" =~ ^[0-9]+$ ]]; then
                echo "  [ERROR] $base 的 N spec 非法: '$nspec' (应为 'auto' 或正整数), 跳过"
                continue
            fi
            echo "        N = $nspec (固定值, 文件无 n 列)"

            python "$MUNGE" \
                --sumstats "$raw" \
                --N "$nspec" \
                --merge-alleles "$MERGE_ALLELES" \
                --out "$out_prefix"
        fi
    fi

    if [ ! -f "$sumstats" ]; then
        echo "  [ERROR] munge 失败, 未生成 $sumstats"
        continue
    fi

    NEW_FILES+=("$sumstats")
    NEW_NAMES+=("$name")
done

NNEW=${#NEW_FILES[@]}
if [ "$NNEW" -eq 0 ]; then
    echo "[FATAL] 没有成功生成任何新 sumstats, 退出"
    exit 1
fi

# ==================== 收集所有"已有"文件 ====================
# 新文件在 MUNGE_OUT_DIR, 不会和下面的扫描冲突
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
echo "  新增 target: $NNEW 个"
echo "  已有文件:    $NEXIST 个"
echo "  待补 rg 对数: 新×旧 $((NNEW * NEXIST)) + 新×新 $((NNEW*(NNEW-1)/2))"
echo "============================================================"
echo ""

# ==================== Step 1: 新文件的 h2 ====================
echo "====== Step 1: 新增 target 的 Heritability ======"

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
        echo "  [ERROR] ${src}:${name} — log 未生成"
        echo "${name},${src},NA,NA,NA,NA" >> "$H2_CSV"
        continue
    fi

    lambda_gc=$(grep "Lambda GC:"    "$log" | awk '{print $NF}')
    mean_chi2=$(grep "Mean Chi^2:"   "$log" | awk '{print $NF}')
    intercept=$(grep "Intercept:"    "$log" | awk '{print $2}')
    intercept_se=$(grep "Intercept:" "$log" | sed 's/.*(\(.*\))/\1/')

    echo "${name},${src},${lambda_gc},${mean_chi2},${intercept},${intercept_se}" >> "$H2_CSV"
done

# ==================== Step 2: 补算 rg ====================
echo ""
echo "====== Step 2: 补算 Genetic Correlation ======"

parse_and_append_rg () {
    local log="$1" name_i="$2" name_j="$3"

    if [ ! -f "$log" ] || ! grep -q "Genetic Correlation:" "$log" 2>/dev/null; then
        echo "        [ERROR] 失败, 查看: $log"
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

# ---- (A) 新文件 vs 所有已有文件 ----
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
            echo "        [SKIP] 已完成"
        else
            python "$LDSC" --rg "${NEW_FILES[$a]},${EXIST_FILES[$b]}" \
                --ref-ld-chr "$REF_LD" \
                --w-ld-chr "$W_LD" \
                --out "$out"
        fi
        parse_and_append_rg "$log" "$name_i" "$name_j"
    done
done

# ---- (B) 新文件之间两两 ----
for (( a=0; a<NNEW; a++ )); do
    for (( b=a+1; b<NNEW; b++ )); do
        ((count++))
        name_i="${NEW_NAMES[$a]}"; src_i="T"
        name_j="${NEW_NAMES[$b]}"; src_j="T"

        out="${LOG_DIR}/${src_i}__${name_i}__vs__${src_j}__${name_j}"
        log="${out}.log"
        echo "  [$count/$NPAIR] ${src_i}:${name_i} vs ${src_j}:${name_j}"

        if [ -f "$log" ] && grep -q "Analysis finished" "$log" 2>/dev/null; then
            echo "        [SKIP] 已完成"
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
echo "  补充分析完成!"
echo "  sumstats 输出: $MUNGE_OUT_DIR/"
echo "  Heritability:  $H2_CSV"
echo "  Genetic Corr:  $RG_CSV"
echo "============================================================"