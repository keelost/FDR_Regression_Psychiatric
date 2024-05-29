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
    level.target = c('magma','smultixcan') 
    head.drug = data.frame(level3_codes=c('N03A','N06A','N05A','N04B','N06B'),NAME=c('Anti-epileptics','Antidepressants','Antipsychotics','Dopaminergic agents','psychostimulants'))
    top.drug = c('level3_codes','NAME','p.adj','target','type','level')
    for (x in target){
        for (y in type.target){
            for (z in level.target){
                drug.detail = fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/',x,y,z,'.category-two.csv'))
                drug.exact = left_join(head.drug,drug.detail[,c(1,5)],by='level3_codes')
                drug.exact$target =x
                drug.exact$type = y
                drug.exact$level = z
                top.drug = rbind(top.drug,drug.exact)
            }
        }
    }
    
    top.drug.p = c('level3_codes','NAME','p.adj','target','type','level')
    for (x in target){
        for (y in type.target){
            for (z in level.target){
                drug.detail = fread(paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/p-version.',x,y,z,'.category-two.csv'))
                drug.exact = left_join(head.drug,drug.detail[,c(1,5)],by='level3_codes')
                drug.exact$target =x
                drug.exact$type = y
                drug.exact$level = z
                top.drug.p = rbind(top.drug.p,drug.exact)
            }
        }
    }
    fwrite(as.data.frame(top.drug),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/5Drug.detail.csv'),sep = ',')
    fwrite(as.data.frame(top.drug.p),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/DEA/p-5Drug.detail.csv'),sep = ',')
    
q('no')
