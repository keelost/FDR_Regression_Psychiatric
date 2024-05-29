#!/bin/bash

cd /exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013
###Part1
for i in *.overlap.4magma.txt
do
brain=('Brain_Amygdala'
'Brain_Anterior_cingulate_cortex_BA24'
'Brain_Caudate_basal_ganglia'
'Brain_Cerebellar_Hemisphere'
'Brain_Cerebellum'
)

for T in "${brain[@]}"
do
	echo "$T"
	nohup /exeh_4/jinghong_qiu/Program/MetaXcan/software/MetaXcan.py \
--model_db_path /exeh_4/jinghong_qiu/Program/MetaXcan/software/data/brain.gtex/gtex_v7_$T'_imputed_europeans_tw_0.5_signif.db' \
--covariance /exeh_4/jinghong_qiu/Program/MetaXcan/software/data/brain.gtex/gtex_v7_$T'_imputed_eur_covariances.txt.gz' \
--gwas_file $i \
--snp_column snpid \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--zscore_column z.decor \
--pvalue_column p.decor \
--output_file /exeh_4/jinghong_qiu/SO_Lab/level2/metaxcan_1_all/mdd2013/gtex_v7_$i'_in_'$T.csv \
>/exeh_4/jinghong_qiu/logfile/mdd2013.brain13.txt &
done
done


###Part2
for i in *.overlap.4magma.txt
do
brain=('Brain_Cortex'
'Brain_Frontal_Cortex_BA9'
'Brain_Hippocampus'
'Brain_Hypothalamus'
)

for T in "${brain[@]}"
do
	echo "$T"
	nohup /exeh_4/jinghong_qiu/Program/MetaXcan/software/MetaXcan.py \
--model_db_path /exeh_4/jinghong_qiu/Program/MetaXcan/software/data/brain.gtex/gtex_v7_$T'_imputed_europeans_tw_0.5_signif.db' \
--covariance /exeh_4/jinghong_qiu/Program/MetaXcan/software/data/brain.gtex/gtex_v7_$T'_imputed_eur_covariances.txt.gz' \
--gwas_file $i \
--snp_column snpid \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--zscore_column z.decor \
--pvalue_column p.decor \
--output_file /exeh_4/jinghong_qiu/SO_Lab/level2/metaxcan_1_all/mdd2013/gtex_v7_$i'_in_'$T.csv \
>/exeh_4/jinghong_qiu/logfile/mdd2013.brain13.txt &
done
done

###Part3
for i in *.overlap.4magma.txt
do
brain=('Brain_Nucleus_accumbens_basal_ganglia'
'Brain_Putamen_basal_ganglia'
'Brain_Spinal_cord_cervical_c-1'
'Brain_Substantia_nigra'
)

for T in "${brain[@]}"
do
	echo "$T"
	nohup /exeh_4/jinghong_qiu/Program/MetaXcan/software/MetaXcan.py \
--model_db_path /exeh_4/jinghong_qiu/Program/MetaXcan/software/data/brain.gtex/gtex_v7_$T'_imputed_europeans_tw_0.5_signif.db' \
--covariance /exeh_4/jinghong_qiu/Program/MetaXcan/software/data/brain.gtex/gtex_v7_$T'_imputed_eur_covariances.txt.gz' \
--gwas_file $i \
--snp_column snpid \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--zscore_column z.decor \
--pvalue_column p.decor \
--output_file /exeh_4/jinghong_qiu/SO_Lab/level2/metaxcan_1_all/mdd2013/gtex_v7_$i'_in_'$T.csv \
>/exeh_4/jinghong_qiu/logfile/mdd2013.brain13.txt &
done
done
