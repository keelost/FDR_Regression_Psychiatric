#************************************#
#-----Get detail of 2 kinds of sample t test-----#
#************************************#
#R-4.0.2
#2020/11/21
#Written by JinghongQIU
#----load packages and set dir----#
	rm(list = ls())
	library(data.table);library(stringr);library(dplyr)
    setwd('/exeh_4/jinghong_qiu/SO_Lab/DEA/')

#read magma results for FDRreg results (zscore-->pvalue-->magma)and rename Pvalue
    target = c('wal.scz2018','bd2018','mdd2019','scz2012','scz2014','scz.eas2019','bd2012','mddco','mdd2013')
    level.target = c('magma','smultixcan')
    
    for (x in target){
        for (y in level.target){
            drug.detail = fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,'qval',y,'.category-two.csv'))
            colnames(drug.detail) = gsub('p.adj','p.adj.two',colnames(drug.detail))
            drug.detail$p.adj.one =999
            p.adj.one = p.adjust(drug.detail$pval.oneSamp.t[which(!drug.detail$pval.oneSamp.t==999)],method='fdr')
            drug.detail$p.adj.one[1:length(p.adj.one)] <- p.adj.one
            drug.detail = drug.detail[,c(1,2,3,6,4,5)]
            drug.detail$target = x
            drug.detail.after = fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,'fdr',y,'.category-two.csv'))
            
            drug.detail.after$p.adj.one =999
            p.adj.one = p.adjust(drug.detail.after$pval.oneSamp.t[which(!drug.detail.after$pval.oneSamp.t==999)],method='fdr')
            drug.detail.after$p.adj.one[1:length(p.adj.one)] <- p.adj.one
            drug.detail.after = drug.detail.after[,c(1,2,3,6,4,5)]
            colnames(drug.detail.after) = c('codes','name','one.t','p.one','two.t','p.two')
            drug.detail = cbind(drug.detail,drug.detail.after)
            fwrite(as.data.frame(drug.detail),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,y,'.cate.all.csv'),sep = ',')
        }
    }

    
q('no')
