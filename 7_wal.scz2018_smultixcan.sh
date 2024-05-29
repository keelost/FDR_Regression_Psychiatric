
b='overlap.4magma.txt_in_(.*).csv'
phenotypes=('ad' 'adhd2019' 'alco2018' 'anx' 'asd2019' 'bd2018' 'cannabis' 'chrono' 'dx.dep' 'dx.life.threat.ill' 'ed2019' 'ever.csh' 'ever.sh' 'ever.st' 'inte' 'life.threat.accident' 'mdd2019' 'mddco' 'neuroti2018' 'no.dep' 'ocd2018' 'pd' 'physical.crime' 'psycho.experi' 'ptsd2019' 'reac' 'sleep.dura' 'smoking' 'st' 'st.dep' 'tense' 'wit.death' 'wal.scz2018')

for T in ${phenotypes[*]}
do
nohup /exeh_3/rstao/programs/MetaXcan/software/SMulTiXcan.py \
--models_folder /exeh_3/rstao/programs/MetaXcan/software/data/brain.gtex/ \
--models_name_pattern 'gtex_v7_(.*)_imputed_europeans_tw_0.5_signif.db' \
--snp_covariance /exeh_3/rstao/programs/multixcan/data/models_all/snp_covariance_v7.txt.gz  \
--metaxcan_folder /exeh_4/jinghong_qiu/SO_Lab/level2/metaxcan_1_all/wal.scz2018/ \
--metaxcan_filter  gtex_v7_${T}.${b} \
--metaxcan_file_name_parse_pattern 'gtex_v7_(.*).overlap.4magma.txt_in_(.*).csv' \
--gwas_file /exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/wal.scz2018/${T}.overlap.txt \
--snp_column snpid \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--zscore_column z \
--pvalue_column pval \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output /exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/${T}.allbrain.txt \
>/exeh_4/jinghong_qiu/logfile/wal.scz2018/${T}.brain12.txt 2>&1 & 
done
