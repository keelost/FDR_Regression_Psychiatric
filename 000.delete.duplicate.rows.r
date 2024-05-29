    rm(list = ls())
	library(data.table);library(stringr);library(Hotelling);library(dplyr);library(mppa);library(ICSNP)
    setwd('/exeh_4/jinghong_qiu/SO_Lab/DEA/')
    load("/exeh_4/jinghong_qiu/SO_Lab/mat.drug_DSigDB.Rdata")
    
    ls()
    test = mat.drug[,-1]
    test2 = distinct(test)
    dup.gene = setdiff(rownames(test),rownames(test2))
    length(dup.gene)
    