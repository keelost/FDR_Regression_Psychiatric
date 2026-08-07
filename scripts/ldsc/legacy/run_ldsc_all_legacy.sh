#!/bin/bash

BASE_DIR="/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all"
LDSC_DIR="/exeh_4/jinghong_qiu/ldsc"
REF_LD="/exeh_4/jinghong_qiu/ldsc/eur_w_ld_chr/"
W_LD="/exeh_4/jinghong_qiu/ldsc/eur_w_ld_chr/"
OUT_DIR="/exeh_4/jinghong_qiu/SO_Lab/15.fdrreg.update/01.result/05.locfdr/ldsc_results"

mkdir -p ${OUT_DIR}

while read TRAIT NCASE NCON; do

    FILE="${BASE_DIR}/${TRAIT}/${TRAIT}.overlap.txt"

    if [ -f "$FILE" ]; then

        echo "Processing ${TRAIT}"

        # Step 1: normalize column names.
        awk 'NR==1{
            print "SNP CHR BP A1 A2 Z INFO P"
            next
        }
        {print $1,$2,$3,$4,$5,$6,$7,$8}' OFS=" " $FILE \
        > ${OUT_DIR}/${TRAIT}.formatted.txt

        # Step 2: munge
        python ${LDSC_DIR}/munge_sumstats.py \
            --sumstats ${OUT_DIR}/${TRAIT}.formatted.txt \
            --N-cas ${NCASE} \
            --N-con ${NCON} \
            --info-min 0.9 \
            --out ${OUT_DIR}/${TRAIT} \
            --merge-alleles ${LDSC_DIR}/w_hm3.snplist

        # Step 3: LDSC
        python ${LDSC_DIR}/ldsc.py \
            --h2 ${OUT_DIR}/${TRAIT}.sumstats.gz \
            --ref-ld-chr ${REF_LD} \
            --w-ld-chr ${W_LD} \
            --out ${OUT_DIR}/${TRAIT}

    else
        echo "File not found: ${FILE}"
    fi

done < /exeh_4/jinghong_qiu/SO_Lab/15.fdrreg.update/01.result/02.mtag/sample_size.txt
