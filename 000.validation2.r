#********************#
#-----Validation-----#
#********************#
#R-4.0.2
#2020/10/15
#Written by JinghongQIU
#----load packages and set dir----#
	rm(list = ls())
	library(data.table);
    setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/')

#----load data for caculation----#  
    head.name=c('type','fdr.sig','fdrr.sig','fold.change','overlap','sig.big','pval')
#_____SNP
#BD
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/bd2018/bd2018.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/bd2012/bd2012.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$snpid,small_data$snpid)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result

#SCZ:SCZ2014
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2014/scz2014.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2012/scz2012.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$snpid,small_data$snpid)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result

#SCZ:WAL.SCZ2018
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/wal.scz2018/wal.scz2018.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2012/scz2012.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$snpid,small_data$snpid)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result

#MDD
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/mdd2019/mdd2019.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/mdd2013/mdd2013.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$snpid[which(big_data$qval<0.05)],small_data$snpid[which(small_data$FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$snpid,small_data$snpid)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result

#_____MAGMA
#BD
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result/bd2018.gene.bio.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2012/Result/bd2012.gene.bio.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ID,small_data$ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    result

#SCZ: SCZ2014
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2014/Result/scz2014.gene.bio.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2012/Result/scz2012.gene.bio.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ID,small_data$ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    result

#SCZ: WAL.SCZ2018
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/wal.scz2018/Result/wal.scz2018.gene.bio.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2012/Result/scz2012.gene.bio.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ID,small_data$ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    result

#MDD
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2019/Result/mdd2019.gene.bio.fdrreg.txt')
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013/Result/mdd2013.gene.bio.fdrreg.txt')
    #get gene number
    fdr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ID[which(big_data$qval<0.05)],small_data$ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ID,small_data$ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    result

#_____PrediXcan
    #set brain regions
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
    #BD
    table.result = c('type','fdr.sig','fdrr.sig','fold.change','overlap','sig.big','pval','region')
    for (region in brain){
    #load data
    big_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2018/Result/bd2018.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    small_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2012/Result/bd2012.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    no.bio = c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value,region)
    bio = c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio,region)
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        lasso = c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso,region)
    } else {
        lasso = c('lasso','NA','NA','NA','NA','NA','NA',region)
    }
    table.result = data.frame(table.result,no.bio, bio,lasso)
    }
    fwrite(as.data.frame(t(table.result)),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/Validation-meta-bd.csv',sep=',')

#SCZ:SCZ2014
    table.result = c('type','fdr.sig','fdrr.sig','fold.change','overlap','sig.big','pval','region')
    for (region in brain){
    #load data
    big_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2014/Result/scz2014.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    small_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/Result/scz2012.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    no.bio = c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value,region)
    bio = c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio,region)
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        lasso = c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso,region)
    } else {
        lasso = c('lasso','NA','NA','NA','NA','NA','NA',region)
    }
    table.result = data.frame(table.result,no.bio, bio,lasso)
    }
    fwrite(as.data.frame(t(table.result)),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/Validation-meta-scz2014.csv',sep=',')

#SCZ:WAL.SCZ2018
    table.result = c('type','fdr.sig','fdrr.sig','fold.change','overlap','sig.big','pval','region')
    for (region in brain){
    #load data
    big_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/Result/wal.scz2018.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    small_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/Result/scz2012.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    no.bio = c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value,region)
    bio = c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio,region)
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        lasso = c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso,region)
    } else {
        lasso = c('lasso','NA','NA','NA','NA','NA','NA',region)
    }
    table.result = data.frame(table.result,no.bio, bio,lasso)
    }
    fwrite(as.data.frame(t(table.result)),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/Validation-meta-wal.scz2018.csv',sep=',')


#MDD
    table.result = c('type','fdr.sig','fdrr.sig','fold.change','overlap','sig.big','pval','region')
    for (region in brain){
    #load data
    big_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2019/Result/mdd2019.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    small_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2013/Result/mdd2013.finalresult.',region,'.test_SetSeed.txt'),sep=' ',header = T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    no.bio = c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value,region)
    bio = c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio,region)
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        lasso = c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso,region)
    } else {
        lasso = c('lasso','NA','NA','NA','NA','NA','NA',region)
    }
    table.result = data.frame(table.result,no.bio, bio,lasso)
    }
    fwrite(as.data.frame(t(table.result)),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/Validation-meta-scz2014.csv',sep=',')

#_____SMULTIXCAN
#BD
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2018/bd2018.finalresult.txt',fill=T)
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2012/bd2012.finalresult.txt',fill=T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        result = cbind(result,c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso))
    } else {
        result = cbind(result,c('lasso','NA','NA','NA','NA','NA','NA'))
    }
    result

#SCZ: SCZ2014
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2014/scz2014.finalresult.txt',fill=T)
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/scz2012.finalresult.txt',fill=T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        result = cbind(result,c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso))
    } else {
        result = cbind(result,c('lasso','NA','NA','NA','NA','NA','NA'))
    }
    result

#SCZ: WAL.SCZ2018
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/wal.scz2018.finalresult.txt',fill=T)
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/scz2012.finalresult.txt',fill=T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        result = cbind(result,c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso))
    } else {
        result = cbind(result,c('lasso','NA','NA','NA','NA','NA','NA'))
    }
    result


#MDD
    #load data
    big_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2019/mdd2019.finalresult.txt',fill=T)
    small_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2013/mdd2013.finalresult.txt',fill=T)
    #get gene number
    fdr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$qval<0.05)])
    fdrr.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$FDR.the<0.05)])
    fdrr.bio.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the<0.05)])
    #caculation:fold change
    fold_change = length(fdrr.sig_gene)/length(fdr.sig_gene)
    fold_change.bio = length(fdrr.bio.sig_gene)/length(fdr.sig_gene)
    overlap = intersect(big_data$ENSEMBL_GENE_ID,small_data$ENSEMBL_GENE_ID)
    sig.big.data = length(which(big_data$qval<0.05))
    matri.data = matrix(c(length(fdr.sig_gene),length(fdrr.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.sig_gene)),nrow=2)
    if (all(matri.data>5)){
        p.cacu = prop.test(x = c(length(fdr.sig_gene),length(fdrr.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu = fisher.test(matri.data,alternative = "less")
    }
    p.value = p.cacu$p.value
    matri.data.bio = matrix(c(length(fdr.sig_gene),length(fdrr.bio.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.sig_gene)),nrow=2)
    if (all(matri.data.bio>5)){
        p.cacu.bio = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
    } else {
        p.cacu.bio = fisher.test(matri.data.bio,alternative = "less")
    }
    p.value.bio = p.cacu.bio$p.value
    #print result
    result = cbind(head.name,c('no',length(fdr.sig_gene),length(fdrr.sig_gene),fold_change,length(overlap),sig.big.data,p.value))
    result = cbind(result,c('bio',length(fdr.sig_gene),length(fdrr.bio.sig_gene),fold_change.bio,length(overlap),sig.big.data,p.value.bio))
    if ('bio.FDR.the.lasso' %in% colnames(small_data)) {
        fdrr.bio.lasso.sig_gene = intersect(big_data$ENSEMBL_GENE_ID[which(big_data$qval<0.05)],small_data$ENSEMBL_GENE_ID[which(small_data$bio.FDR.the.lasso<0.05)])
        fold_change.bio.lasso = length(fdrr.bio.lasso.sig_gene)/length(fdr.sig_gene)
        matri.data.bio.lasso = matrix(c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),sig.big.data-length(fdr.sig_gene),sig.big.data-length(fdrr.bio.lasso.sig_gene)),nrow=2)
        if (all(matri.data.bio.lasso>5)){
            p.cacu.bio.lasso = prop.test(x = c(length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene)), n = c(sig.big.data, sig.big.data),alternative = "less")
            } else {
                p.cacu.bio.lasso = fisher.test(matri.data.bio.lasso,alternative = "less")
                }
        p.value.bio.lasso = p.cacu.bio.lasso$p.value
        result = cbind(result,c('lasso',length(fdr.sig_gene),length(fdrr.bio.lasso.sig_gene),fold_change.bio.lasso,length(overlap),sig.big.data,p.value.bio.lasso))
    } else {
        result = cbind(result,c('lasso','NA','NA','NA','NA','NA','NA'))
    }
    result

