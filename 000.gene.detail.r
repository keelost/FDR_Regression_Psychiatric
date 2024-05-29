#*****************************#
#-----Detail of sig.genes-----#
#*****************************#
#R-4.0.2
#2020/11/02
#Written by JinghongQIU
#----load packages and set dir----#
	rm(list = ls())
	library(data.table);library(stringr);library(dplyr)
    setwd('/exeh_4/jinghong_qiu/SO_Lab/gene.detail/')
    lib.data = fread('/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-entrez.csv')
#read magma results for FDRreg results (zscore-->pvalue-->magma)and rename Pvalue
    target = c('wal.scz2018','bd2018','mdd2019','bd2012','mddco','mdd2013','scz2012','scz2014','scz.eas2019')
    level.target = c('magma','smultixcan')
    for (x in target){
        for (y in level.target){
            cat(x,y,'is working now\n')
            if ('magma' %in% y){
                trait_gene=fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/',x,'/Result/',x,'.gene.bio.fdrreg.txt'))
                trait_gene = left_join(trait_gene,lib.data[,c(1,2)],by='ID')
            } else {
                trait_gene=read.table(file=paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/',x,'/',x,'.finalresult.txt'),sep = ' ',header=T)                
            }
            if (length(which(trait_gene$bio.FDR.the<0.05))>0){
                sig.target = trait_gene[which(trait_gene$bio.FDR.the<0.05),]
                fwrite(as.data.frame(sig.target),file=paste0(x,'-',y,'.csv'),sep=',')
            } else {
                cat('no sig genes\n')
            }
            cat(x,y,'is finished\n')
        }
        #sig.gene in 13 brain regions
        brain=c('Brain_Amygdala',
        'Brain_Anterior_cingulate_cortex_BA24',
        'Brain_Caudate_basal_ganglia',
        'Brain_Cerebellar_Hemisphere',
        'Brain_Cerebellum',
        'Brain_Cortex',
        'Brain_Frontal_Cortex_BA9',
        'Brain_Hippocampus',
        'Brain_Hypothalamus',
        'Brain_Nucleus_accumbens_basal_ganglia',
        'Brain_Putamen_basal_ganglia',
        'Brain_Spinal_cord_cervical_c-1',
        'Brain_Substantia_nigra'
        )
        table.name = c("ENSEMBL_GENE_ID","ENSEMBL_GENE_ID_name","zscore","effect_size","pvalue","var_g","pred_perf_r2","pred_perf_pval","pred_perf_qval","n_snps_used","n_snps_in_cov","n_snps_in_model","qval","FDR.the","FDR.emp","bio.FDR.the","bio.FDR.emp",'region')
        for (i in brain){        
            cat(x,i,'is working\n')
            trait_gene=read.table(file=paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/',x,'/Result/',x,'.finalresult.',i,'.test_SetSeed.txt'),sep = ' ',header=T)
            trait_gene = trait_gene[,c(1:17)]
            trait_gene$region = i
            if (length(which(trait_gene$bio.FDR.the<0.05))>0) {
                sig.target = trait_gene[which(trait_gene$bio.FDR.the<0.05),]
                table.name = rbind(table.name,sig.target)
                cat(x,i,'is finished\n')
            } else {
                cat(x,i,'is not sig\n')
            }
            
        }
        fwrite(as.data.frame(table.name),file=paste0(x,'-meta','.csv'),sep=',')
    }
    
q('no')
