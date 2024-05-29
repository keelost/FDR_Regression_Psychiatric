
b='overlap.4magma.txt_in_(.*).csv'
phenotypes=('alco2018' 'asd2019' 'bd2018' 'cannabis' 'chrono' 'ed2019' 'edu' 'ever.st' 'insom2019' 'inte' 'life.threat.accident' 'mdd2019' 'mddco' 'neuroti2018' 'no.dep' 'ocd2018' 'physical.crime' 'ptsd2019' 'reac' 'sleep.dura')

for T in ${phenotypes[*]}
do
nohup /exeh_3/rstao/programs/MetaXcan/software/SMulTiXcan.py \
--models_folder /exeh_3/rstao/programs/MetaXcan/software/data/brain.gtex/ \
--models_name_pattern 'gtex_v7_(.*)_imputed_europeans_tw_0.5_signif.db' \
--snp_covariance /exeh_3/rstao/programs/multixcan/data/models_all/snp_covariance_v7.txt.gz  \
--metaxcan_folder /exeh_4/jinghong_qiu/SO_Lab/level2/metaxcan_1_all/scz.eas2019/ \
--metaxcan_filter  gtex_v7_${T}.${b} \
--metaxcan_file_name_parse_pattern 'gtex_v7_(.*).overlap.4magma.txt_in_(.*).csv' \
--gwas_file /exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz.eas2019/${T}.overlap.txt \
--snp_column snpid \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--zscore_column z \
--pvalue_column pval \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output /exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz.eas2019/${T}.allbrain.txt \
>/exeh_4/jinghong_qiu/logfile/scz.eas2019/${T}.brain12.txt 2>&1 & 
done
