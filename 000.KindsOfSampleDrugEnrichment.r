#************************************#
#-----Get detail info of 5 drugs-----#
#************************************#
#R-4.0.2
#2020/11/17
#Written by JinghongQIU
#----load packages and set dir----#
	rm(list = ls())
	library(data.table);library(stringr);library(dplyr)
    setwd('/exeh_4/jinghong_qiu/SO_Lab/DEA/')

#read magma results for FDRreg results (zscore-->pvalue-->magma)and rename Pvalue
    target = c('wal.scz2018','bd2018','mdd2019','scz2012','scz2014','scz.eas2019','bd2012','mddco','mdd2013')
    type.target = c('qval','fdr')
    head.drug = data.frame(level3_codes=c('N03A','N06A','N05A','N04B','N06B','N05B','N05C'),NAME=c('Anti-epileptics','Antidepressants','Antipsychotics','Dopaminergic agents','Psychostimulants','anxiolytic','hypnotics'))
    top.drug.one = c('level3_codes','t','p.adj','target','type')
    top.drug.two = top.drug.one
    for (x in target){
        for (y in type.target){
            drug.detail = fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,y,'magma.category-two.csv'))
            drug.detail$p.adj.one =999
            p.adj.one = p.adjust(drug.detail$pval.oneSamp.t[which(!drug.detail$pval.oneSamp.t==999)],method='fdr')
            drug.detail$p.adj.one[1:length(p.adj.one)] <- p.adj.one
            drug.detail = drug.detail[,c(1,2,3,6,4,5)]
            drug.exact = left_join(head.drug,drug.detail,by='level3_codes')
            drug.exact$target =x
            drug.exact$type = y
            drug.exact.one = drug.exact[,c(1,4,5,8,9)]
            drug.exact.two = drug.exact[,c(1,6,7,8,9)]
            top.drug.one = rbind(top.drug.one,drug.exact.one)
            top.drug.two = rbind(top.drug.two,drug.exact.two)
        }
    }
    fwrite(as.data.frame(top.drug.one),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/7Drugone.detail.csv'),sep = ',')
    fwrite(as.data.frame(top.drug.two),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/7Drugtwo.detail.csv'),sep = ',')
    
q('no')
