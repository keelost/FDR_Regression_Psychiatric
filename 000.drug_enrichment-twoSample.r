#*************************#
#-----Drug Enrichment-----#
#*************************#
#R-4.0.2
#2020/11/13
#Written by JinghongQIU
#----load packages and set dir----#
	rm(list = ls())
	library(data.table);library(stringr);library(Hotelling);library(dplyr);library(mppa);library(ICSNP)
    setwd('/exeh_4/jinghong_qiu/SO_Lab/DEA/')
    source('/exeh_3/rstao/dr.so/002.meta/003.code/fdr.reg.ptsd/drug.enrichment.analysis_function.r')

#read magma results for FDRreg results (zscore-->pvalue-->magma)and rename Pvalue
    lib.data = fread('/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-entrez.csv')
    target = c('wal.scz2018','bd2018','mdd2019','scz2012','scz2014','scz.eas2019','bd2012','mddco','mdd2013')
    type.target = c('qval','fdr')
    level.target = c('magma','smultixcan') 
    for (x in target){
        for (z in level.target){
            if ('magma' %in% z){
                trait_gene=fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/',x,'/Result/',x,'.gene.bio.fdrreg.txt'))
                #check if gene.symbol exist
                if (!'GeneName' %in% colnames(trait_gene)){trait_gene = left_join(trait_gene,lib.data[,c(1,2)],by='ID')}
                colnames(trait_gene)[colnames(trait_gene)=='GeneName'] <- "Gene"
            } else {
                trait_gene=read.table(file=paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/',x,'/',x,'.finalresult.txt'),sep = ' ',header=T)
                #check if gene.symbol exist
                colnames(trait_gene)[colnames(trait_gene)=='ENSEMBL_GENE_ID_name'] <- "Gene"
            }
    
    for (y in type.target){    
        cat(x,y,z,'is working now.\n')
        if ('qval' %in% y){    
            colnames(trait_gene)[colnames(trait_gene)=='qval'] <- "Pvalue"} else {        
                colnames(trait_gene)[colnames(trait_gene)=='bio.FDR.the'] <- "Pvalue"    
            }
#load drug database
    load("/exeh_4/jinghong_qiu/SO_Lab/mat.drug_DSigDB.Rdata")
    drug_enrichment_test <- function(GWAS_GeneBasedTest_file = trait_gene,gene_col = trait_gene$Gene,pvalue_col = trait_gene$Pvalue ) {
        attach(trait_gene)
        ## set pvalue which are equal to 1 to be 0.9999, as qnorm(1) will return Inf
        Pvalue[Pvalue == 1] <- 0.9999
        mat.disease = data.frame(Gene = Gene,zval_trait = qnorm(Pvalue))
        detach(trait_gene)

        intersect.genes = intersect(mat.drug$Gene, mat.disease$Gene)

        mat.disease2  = mat.disease[which(mat.disease$Gene %in% intersect.genes), ]
        mat.disease2  = arrange(mat.disease2, Gene)
        #some z will become -Inf, change them into -5e-08
        mat.disease2$zval_trait[which(mat.disease2$zval_trait==-Inf)] <- -5e-08
        mat.drug2 = mat.drug[which(mat.drug$Gene %in% intersect.genes),]
        mat.drug2  = arrange(mat.drug2, Gene)
        
        #there are some duplicated gene stopping the t-test
        if (!length(mat.disease2$Gene)==length(mat.drug2$Gene)){        
            mat.disease2 = mat.disease2[-which(duplicated(mat.disease2$Gene)),]
            }

        t1 = proc.time()
        no.drugs = ncol(mat.drug2)-1
        t.test.p.one.sided_trait = as.numeric(no.drugs)

        #********Start of loop*********************
        for (i in 1:no.drugs ) { 
            genes.for.drug = mat.drug2[,i+1]  
            no.genes.for.drug = sum(genes.for.drug==1)

        ##note that we only consider a drug for analysis if there are at least 5 associated genes
            if (no.genes.for.drug <5) { 
                t.test.p.one.sided_trait[i]=999; 
                } else {
                t.test.obj = t.test(mat.disease2$zval_trait ~ genes.for.drug, alternative="greater")
                t.test.p.one.sided_trait[i] <- t.test.obj$p.value
            }
        }
        t2=proc.time()
        t2-t1
        #*********************

        res = data.frame(drug= colnames(mat.drug2)[-1],
                        t.test.p.one.sided_trait = t.test.p.one.sided_trait) 
        return(res)				 
        }

#perform drug enrichment analysis for magma results
    DEA=drug_enrichment_test(trait_gene,trait_gene$Gene,trait_gene$Pvalue)
    DEA=DEA[order(DEA$t.test.p.one.sided_trait),]
    DEA$p.adj = 999
    p.adj = p.adjust(DEA$t.test.p.one.sided_trait[which(!DEA$t.test.p.one.sided_trait==999)],method='fdr')
    DEA$p.adj[1:length(p.adj)] <- p.adj
    
#save results of DEA
    save(DEA,file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,y,z,'.enrich.rdata'))
    fwrite(as.data.frame(DEA),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,y,z,'.enrich.csv'),sep = ',')

    load('/exeh_3/rstao/dr.so/002.meta/002.clear.data/ATC_drug_lists_ALL.Rdata') #load all ATC categories
    source('/exeh_3/rstao/dr.so/002.meta/003.code/sa.pgc/005.drug_ATC.category_asFunc.R')

#load results from drug enrichment analysis
#target='hippocampal_vol'
    result_file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,y,z,'.enrich.csv')
    pvalue_col_name="t.test.p.one.sided_trait"
    drug_col_name='drug'
    ATC_enrich_res = ATC_enrichment(result_file,pvalue_col_name,drug_col_name)
    #Attention: it is the twoSamp version
    ATC_enrich_res=ATC_enrich_res[order(ATC_enrich_res$pval.twoSamp.t),]
    ATC_enrich_res$p.adj = 999
    p.adj.atc = p.adjust(ATC_enrich_res$pval.twoSamp.t[which(!ATC_enrich_res$pval.twoSamp.t==999)],method='fdr')
    ATC_enrich_res$p.adj[1:length(p.adj.atc)] <- p.adj.atc
#write.csv(ATC_enrich_res,"ATC_all_enrichment.csv", quote=F, row.names=F)
    fwrite(as.data.frame(ATC_enrich_res),file=paste0(x,y,z,'.category-two.csv'),sep=',')
    cat(x,y,z,'is finished.\n')
    }}}
q('no')
