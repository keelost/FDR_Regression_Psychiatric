#***********************************#
#-----FDR_REG: mdd2019-metaxcan-----#
#***********************************#
#R-4.0.2
#2020/09/07
#Written by JinghongQIU

#----load packages and set dir----#
	rm(list = ls())
	library(FDRreg);library(data.table);library(powerplus);library(stringr);library(dplyr);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/metaxcan_1_all/mdd2019')
#get biological information from online resource and merge with gene list#
#1.DAVID Bioinformatics#
	scz.gene.david=fread(file ='/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-ensembl.csv')
#----Loading data----#
#Set brain regions
	#brain_region = 'Brain_Anterior_cingulate_cortex_BA24'
	#'Brain_Substantia_nigra'
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
	final_sum = c(1:9)
	range_limit = c('<0.5','<0.4','<0.3','<0.2','<0.1','<0.05','<0.04','<0.03','<0.02','<0.01','<0.001','<5e-04','<5e-06','<5e-08')
	assessment.final.sum = c(1:11)
for (brain_region in brain){

#read all predixcan results in cerebellum region
	filenames=list.files(path='.',pattern = brain_region,full.names = T)
	files=lapply(filenames,fread,stringsAsFactor=F)
	target_rank = which(str_detect(filenames,'mdd2019'))
	for (i in 1:length(files)) {
	files[[i]]=files[[i]][complete.cases(files[[i]]$pvalue),]
	files[[i]]=files[[i]][order(files[[i]]$gene),]
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
	target.z=files[[target_rank]]$zscore
	colnames(z.combine)=c('ad','adhd2019','alco2018','antisocial','anx','asd2019','bd2018','cannabis','chrono','college2013','dx.dep','dx.life.threat.ill','ed2019','edu','ever.csh','ever.sh','ever.st','insom2019','inte','life.threat.accident','neuroti2018','no.dep','ocd2018','physical.crime','psycho.experi','ptsd2019','sa.ipsych','scz.eas2019','sleep.dura','st','st.dep','tense','wal.scz2018','wit.death')

#generate geneID and output#
	gene.list=files[[target_rank]]$gene
	gene.list2 = strsplit(gene.list, '[.]')
	listdata = vector()
	for (line in 1:length(gene.list2) ){
	listdata[line] = gene.list2[[line]][1]
	}
	
#perform FDRreg with two different models
	set.seed(200)
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
	set.seed(200)
	fdr.target.theore=FDRreg(target.z,z.combine.bio,nulltype = 'theoretical',method = 'pr')
	fdr.target.empiri=FDRreg(target.z,z.combine.bio,nulltype = 'empirical',method = 'pr')
	
#assessment-fdr.reg-bio
	 features.se=SEfromHessian(fdr.target.theore$model$hessian)
	 features.coef=fdr.target.theore$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment2 = data.frame(colnames(z.combine.bio),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment2) = c('feature','p2','beta2','se2')
	
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
	fdr.bio.the <- c(length(which(target.gene.order$bio.FDR.the<0.5)),
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
	fdr.bio.emp <- c(length(which(target.gene.order$bio.FDR.emp<0.5)),
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
	assessment.final = merge(assessment,assessment2,by='feature',all=T)
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
	 assessment.final = merge(assessment.final,assessment3,by='feature',all=T)
	 set.seed(200)
	 target.gene.order$bio.FDR.the.lasso = lasso.fdr.target.theore$FDR
	 target.gene.order$bio.FDR.emp.lasso = lasso.fdr.target.empiri$FDR
	 fdr.bio.the.lasso <- c(length(which(target.gene.order$bio.FDR.the.lasso<0.5)),
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
	 fdr.bio.emp.lasso <- c(length(which(target.gene.order$bio.FDR.emp.lasso<0.5)),
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
	fdr.bio.the.lasso <- c('NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA')
	fdr.bio.emp.lasso <- c('NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA')
	lasso.fdr.target.theore <- 'NA'
	lasso.fdr.target.empiri <- 'NA'
	assessment.final$p3 = 'NA'
	assessment.final$beta3 = 'NA'
	assessment.final$se3 = 'NA'
}	

#Save final result
	sum_data = data.frame(qval,fdr.the,fdr.bio.the,fdr.bio.the.lasso,fdr.emp,fdr.bio.emp,fdr.bio.emp.lasso)
	sum_data$region = brain_region
	sum_data$range_limit = range_limit
	final_sum = rbind(final_sum,sum_data)
	fwrite(as.data.frame(sum_data),paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/meta_mdd2019_',brain_region,'.test_SetSeed.csv'),sep=',')
	save(fdr.theore,fdr.empiri,fdr.target.theore,fdr.target.empiri,lasso.fdr.target.theore,lasso.fdr.target.empiri, file =paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/meta_mdd2019_',brain_region,'.test_SetSeed.RData'))
	fwrite(as.data.frame(target.gene.order),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2019/Result/mdd2019.finalresult.',brain_region,'.test_SetSeed.txt'),sep=' ')
	assessment.final$region = brain_region
	assessment.final.sum = rbind(assessment.final.sum,assessment.final)
	fwrite(as.data.frame(assessment.final),file=paste0('/exeh_4/jinghong_qiu/SO_Lab/Contribution/meta-mdd2019-',brain_region,'.test_SetSeed.csv'),sep=',')
}
	fwrite(as.data.frame(final_sum),file='/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/meta_mdd2019.all.test_SetSeed200.csv',sep=',')
	fwrite(as.data.frame(assessment.final.sum),file='/exeh_4/jinghong_qiu/SO_Lab/Contribution/meta-mdd2019-all.test_SetSeed200.csv',sep=',')
#exit the software
	q('no') 