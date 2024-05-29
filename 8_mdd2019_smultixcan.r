#*********************************#
#-----FDR_REG: mdd2019-smulti-----#
#*********************************#
#R-4.0.2
#2020/09/18
#Written by JinghongQIU

#----load packages and set dir----#
	rm(list = ls())
	library(FDRreg);library(data.table);library(powerplus);library(stringr);library(dplyr);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2019')
#get biological information from online resource and merge with gene list#
#1.DAVID Bioinformatics#
	scz.gene.david=fread(file ='/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-ensembl.csv')
#----Loading data----#
#read all predixcan results in cerebellum region
	filenames=list.files(path='.',pattern = 'allbrain',full.names = T)
	files=lapply(filenames,fread,stringsAsFactor=F)
	target_rank = which(str_detect(filenames,'mdd2019'))
#delete NA in pvalue column, generate zscore from smultixcan pvalue and order datasets with gene ID
	for (i in 1:length(files)) {
	files[[i]]=files[[i]][complete.cases(files[[i]]$pvalue),]
	files[[i]]=files[[i]][order(files[[i]]$gene),]
	files[[i]]$zscore = qnorm(files[[i]]$pvalue/2)
	}
	other_file = c(2:length(files))
	other_file =other_file[-(target_rank-1)]
#delete NA in pvalue column, generate zscore from smultixcan pvalue and order datasets with gene ID
	if (!identical(files[[target_rank]]$gene,files[[length(files)]]$gene)){
	idx = match(files[[target_rank]]$gene,files[[length(files)]]$gene)
	omit_rank = setdiff(seq(1:tail(idx,n=1)),idx)
	files[[1]] = files[[1]][-omit_rank,]
	for (i in other_file){
	files[[i]] = files[[i]][-omit_rank,]
	}
	}
#extract target zscore and construct variables z.score matrix 
	z.combine = abs(files[[1]]$zscore)
	for (i in other_file){
	z.combine = cbind(z.combine,abs(files[[i]]$zscore))
	}
	colnames(z.combine)=c('ad','adhd2019','alco2018','antisocial','anx','asd2019','bd2018','cannabis','chrono','college2013','dx.dep','dx.life.threat.ill','ed2019','edu','ever.csh','ever.sh','ever.st','insom2019','inte','life.threat.accident','neuroti2018','no.dep','ocd2018','physical.crime','psycho.experi','ptsd2019','sa.ipsych','scz.eas2019','sleep.dura','st','st.dep','tense','wal.scz2018','wit.death')
	set.seed(127)
	target.z=abs(qnorm(files[[target_rank]]$pvalue/2))*sign(rnorm(length(files[[target_rank]]$pvalue)))

#generate geneID and output#
	gene.list=files[[target_rank]]$gene
	gene.list2 = strsplit(gene.list, '[.]')
	listdata = vector()
	for (line in 1:length(gene.list2) ){
	listdata[line] = gene.list2[[line]][1]
	}
	fwrite(as.data.frame(listdata),file = paste0('smultixcan.gene.list.txt'))
	
#perform FDRreg with two different models
	set.seed(127)
	fdr.theore=FDRreg(target.z,z.combine,nulltype = 'theoretical',method = 'pr')
	fdr.empiri=FDRreg(target.z,z.combine,nulltype = 'empirical',method = 'pr')
	
#assessment-fdr-reg
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.theore$model$hessian)
	 features.coef=fdr.theore$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment = data.frame(colnames(z.combine),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment) = c('feature','p1','beta1','se1')
	 
#############################
#---Biological info merge---#
#############################
	target.gene.order = files[[target_rank]]
	target.gene.order$gene = listdata
	colnames(target.gene.order) = gsub("gene","ENSEMBL_GENE_ID",colnames(target.gene.order))
	scz.bio=left_join(target.gene.order[,1:2],scz.gene.david,by='ENSEMBL_GENE_ID')
	scz.bio[is.na(scz.bio)]=0
	
#---------Lasso---------#
	bio.use = as.matrix(scz.bio[,c(3:23)])
	lasso.cv.bio <- cv.glmnet(bio.use,abs(target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
	coef_list.min.bio <- unlist(coef_list.min.bio)[-1] 
	z.combine.lasso.min.bio=bio.use[,coef_list.min.bio != 0]
	
	z.combine.bio=as.matrix(cbind(abs(z.combine),abs(bio.use)))
	set.seed(127)
	fdr.target.theore=FDRreg(target.z,z.combine.bio,nulltype = 'theoretical',method = 'pr')
	fdr.target.empiri=FDRreg(target.z,z.combine.bio,nulltype = 'empirical',method = 'pr')
	
#assessment-fdr.reg-bio
	 features.se=SEfromHessian(fdr.target.theore$model$hessian)
	 features.coef=fdr.target.theore$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment2 = data.frame(colnames(z.combine.bio),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment2) = c('feature','p2','beta2','se2')
	 assessment_sum = merge(assessment,assessment2,by='feature',all=T)

#completion of analysis
	target.gene.order$qval = p.adjust(files[[target_rank]]$pvalue,method='fdr',n=length(files[[target_rank]]$pvalue))
	target.gene.order$FDR.the = fdr.theore$FDR
	target.gene.order$FDR.emp = fdr.empiri$FDR
	target.gene.order$bio.FDR.the = fdr.target.theore$FDR
	target.gene.order$bio.FDR.emp = fdr.target.empiri$FDR

#Show the result <- c()
	qval <- c(length(which(target.gene.order$qval<0.5)),
	length(which(target.gene.order$qval<0.4)),
	length(which(target.gene.order$qval<0.3)),
	length(which(target.gene.order$qval<0.2)),
	length(which(target.gene.order$qval<0.1)),
	length(which(target.gene.order$qval<0.05)),
	length(which(target.gene.order$qval<0.04)),
	length(which(target.gene.order$qval<0.03)),
	length(which(target.gene.order$qval<0.02)),
	length(which(target.gene.order$qval<0.01)),
	length(which(target.gene.order$qval<0.001)),
	length(which(target.gene.order$qval<5e-04)),
	length(which(target.gene.order$qval<5e-06)),
	length(which(target.gene.order$qval<5e-08)))
	fdr.the <- c(length(which(target.gene.order$FDR.the<0.5)),
	length(which(target.gene.order$FDR.the<0.4)),
	length(which(target.gene.order$FDR.the<0.3)),
	length(which(target.gene.order$FDR.the<0.2)),
	length(which(target.gene.order$FDR.the<0.1)),
	length(which(target.gene.order$FDR.the<0.05)),
	length(which(target.gene.order$FDR.the<0.04)),
	length(which(target.gene.order$FDR.the<0.03)),
	length(which(target.gene.order$FDR.the<0.02)),
	length(which(target.gene.order$FDR.the<0.01)),
	length(which(target.gene.order$FDR.the<0.001)),
	length(which(target.gene.order$FDR.the<5e-04)),
	length(which(target.gene.order$FDR.the<5e-06)),
	length(which(target.gene.order$FDR.the<5e-08)))
	fdr.emp <- c(length(which(target.gene.order$FDR.emp<0.5)),
	length(which(target.gene.order$FDR.emp<0.4)),
	length(which(target.gene.order$FDR.emp<0.3)),
	length(which(target.gene.order$FDR.emp<0.2)),
	length(which(target.gene.order$FDR.emp<0.1)),
	length(which(target.gene.order$FDR.emp<0.05)),
	length(which(target.gene.order$FDR.emp<0.04)),
	length(which(target.gene.order$FDR.emp<0.03)),
	length(which(target.gene.order$FDR.emp<0.02)),
	length(which(target.gene.order$FDR.emp<0.01)),
	length(which(target.gene.order$FDR.emp<0.001)),
	length(which(target.gene.order$FDR.emp<5e-04)),
	length(which(target.gene.order$FDR.emp<5e-06)),
	length(which(target.gene.order$FDR.emp<5e-08)))
	bio.fdr.the <- c(length(which(target.gene.order$bio.FDR.the<0.5)),
	length(which(target.gene.order$bio.FDR.the<0.4)),
	length(which(target.gene.order$bio.FDR.the<0.3)),
	length(which(target.gene.order$bio.FDR.the<0.2)),
	length(which(target.gene.order$bio.FDR.the<0.1)),
	length(which(target.gene.order$bio.FDR.the<0.05)),
	length(which(target.gene.order$bio.FDR.the<0.04)),
	length(which(target.gene.order$bio.FDR.the<0.03)),
	length(which(target.gene.order$bio.FDR.the<0.02)),
	length(which(target.gene.order$bio.FDR.the<0.01)),
	length(which(target.gene.order$bio.FDR.the<0.001)),
	length(which(target.gene.order$bio.FDR.the<5e-04)),
	length(which(target.gene.order$bio.FDR.the<5e-06)),
	length(which(target.gene.order$bio.FDR.the<5e-08)))
	bio.fdr.emp <- c(length(which(target.gene.order$bio.FDR.emp<0.5)),
	length(which(target.gene.order$bio.FDR.emp<0.4)),
	length(which(target.gene.order$bio.FDR.emp<0.3)),
	length(which(target.gene.order$bio.FDR.emp<0.2)),
	length(which(target.gene.order$bio.FDR.emp<0.1)),
	length(which(target.gene.order$bio.FDR.emp<0.05)),
	length(which(target.gene.order$bio.FDR.emp<0.04)),
	length(which(target.gene.order$bio.FDR.emp<0.03)),
	length(which(target.gene.order$bio.FDR.emp<0.02)),
	length(which(target.gene.order$bio.FDR.emp<0.01)),
	length(which(target.gene.order$bio.FDR.emp<0.001)),
	length(which(target.gene.order$bio.FDR.emp<5e-04)),
	length(which(target.gene.order$bio.FDR.emp<5e-06)),
	length(which(target.gene.order$bio.FDR.emp<5e-08)))
	
	cat('----------------------- \n All gene number is: ',length(gene.list),sep = "\n")

#If lasso can be applied: 
	if (!is.null(intersect(colnames(bio.use),colnames(z.combine.lasso.min.bio)))){
	z.combine.bio.lasso=as.matrix(cbind(abs(z.combine),abs(z.combine.lasso.min.bio)))
	lasso.fdr.target.theore=FDRreg(target.z,z.combine.bio.lasso,nulltype = 'theoretical',method = 'pr')
	lasso.fdr.target.empiri=FDRreg(target.z,z.combine.bio.lasso,nulltype = 'empirical',method = 'pr')
#assessment-fdr.reg-lasso
	 features.se=SEfromHessian(lasso.fdr.target.theore$model$hessian)
	 features.coef=lasso.fdr.target.theore$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment3 = data.frame(colnames(z.combine.bio.lasso),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment3) = c('feature','p3','beta3','se3')
	 assessment_sum = merge(assessment_sum,assessment3,by='feature',all=T)
	 set.seed(137)
	 target.gene.order$bio.FDR.the.lasso = lasso.fdr.target.theore$FDR
	 target.gene.order$bio.FDR.emp.lasso = lasso.fdr.target.empiri$FDR
	 bio.fdr.the.lasso <- c(length(which(target.gene.order$bio.FDR.the.lasso<0.5)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.4)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.3)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.2)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.1)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.05)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.04)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.03)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.02)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.01)),
	length(which(target.gene.order$bio.FDR.the.lasso<0.001)),
	length(which(target.gene.order$bio.FDR.the.lasso<5e-04)),
	length(which(target.gene.order$bio.FDR.the.lasso<5e-06)),
	length(which(target.gene.order$bio.FDR.the.lasso<5e-08)))
	bio.fdr.emp.lasso <- c(length(which(target.gene.order$bio.FDR.emp.lasso<0.5)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.4)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.3)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.2)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.1)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.05)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.04)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.03)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.02)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.01)),
	length(which(target.gene.order$bio.FDR.emp.lasso<0.001)),
	length(which(target.gene.order$bio.FDR.emp.lasso<5e-04)),
	length(which(target.gene.order$bio.FDR.emp.lasso<5e-06)),
	length(which(target.gene.order$bio.FDR.emp.lasso<5e-08)))
	} else {
	cat('-----\nNo need to lasso\n-----\n')
	assessment_sum$p3 = 'NA'
	assessment_sum$beta3 = 'NA'
	assessment_sum$se3 = 'NA'
	bio.fdr.the.lasso <- c('NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA')
	bio.fdr.emp.lasso <- c('NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA')
	fdr.target.theore.bio.lasso <-'NA'
	fdr.target.empiri.bio.lasso <-'NA'
}
#Save final result
	fwrite(as.data.frame(target.gene.order),file=paste0('/exeh_4/jinghong_qiu/mdd2019.finalresult.test3.txt'),sep=' ')
	fwrite(as.data.frame(assessment_sum),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/smulti-mdd2019.csv',sep=',')
	sum_data = data.frame(qval,fdr.the,fdr.emp,bio.fdr.the,bio.fdr.emp,bio.fdr.the.lasso,bio.fdr.emp.lasso)
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/smulti_mdd2019.csv',sep=',')	
	save(fdr.theore,fdr.empiri,fdr.target.theore,fdr.target.empiri,fdr.target.theore.bio.lasso,fdr.target.empiri.bio.lasso,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/smulti_mdd2019.RData")
#exit the software
	q('no') 