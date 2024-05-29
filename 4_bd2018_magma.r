#*******************************#
#-----FDR_REG: bd2018-magma-----#
#*******************************#
#R-4.0.2
#2020/09/18
#Written by JinghongQIU

#----load packages and set dir----#
	rm(list = ls())
	library(FDRreg);library(data.table);library(powerplus);library(stringr);library(dplyr);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result')

#----Loading data----#

#Target
	target.gene=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result/bd2018.genes.out')
	target.gene.order=target.gene[order(target.gene$GENE),]
#Variables
	listnames = c('ad','adhd2019','alco2018','asd2019','cannabis','college2013','dx.dep','ed2019','edu','ever.csh','ever.sh','ever.st','life.threat.accident','mdd2019','mddco','neuroti2018','no.dep','ocd2018','physical.crime','psycho.experi','ptsd2019','reac','scz.eas2019','sleep.dura','tense','wal.scz2018','wit.death')
	files=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result/',listnames,'.genes.out'),fread,stringsAsFactor=F)
#data preproccess
	fileout.order = files[[1]][order(files[[1]]$GENE),]
	z.combine = abs(qnorm(fileout.order$P/2))
	for (i in 2:length(files)){
	geneout = files[[i]]
	geneout.order=geneout[order(geneout$GENE),]
	z.combine = cbind(z.combine,abs(qnorm(geneout.order$P/2)))
	}
	colnames(z.combine) = listnames

#Simulate random
	set.seed(25)

#extract target z.stat
	target.z=abs(qnorm(target.gene.order$P/2))*sign(rnorm(length(target.gene.order$P)))

#generate geneID and output for biological annotation#
	gene.list=target.gene.order$GENE
	fwrite(as.data.frame(gene.list),'bd2018.gene.list.txt')

#perform FDR regression of 2 nulltypes
	set.seed(25)
	fdr.target.theore=FDRreg(target.z,z.combine,nulltype = 'theoretical',method = 'pr')
	fdr.target.empiri=FDRreg(target.z,z.combine,nulltype = 'empirical',method = 'pr')

#assessment-fdr-reg
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.target.theore$model$hessian)
	 features.coef=fdr.target.theore$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment = data.frame(colnames(z.combine),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment) = c('feature','p1','beta1','se1')

#######################
#--NEW PART-----------#
#######################

#get biological information from online resource and merge with gene list#
#1.DAVID Bioinformatics#
	scz.gene.david=fread('/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-entrez.csv')
	colnames(target.gene.order) = gsub("GENE","ID",colnames(target.gene.order))
	scz.bio=left_join(target.gene.order[,1:3],scz.gene.david,by='ID')
	scz.bio[is.na(scz.bio)]=0
	scz.bio.use = scz.bio[,c(5:25)]
	z.combine.bio=as.matrix(cbind(abs(z.combine),abs(scz.bio[,c(5:25)])))
	set.seed(25)
	fdr.target.theore.bio=FDRreg(target.z,z.combine.bio,nulltype = 'theoretical',method = 'pr')
	fdr.target.empiri.bio=FDRreg(target.z,z.combine.bio,nulltype = 'empirical',method = 'pr')
	
#---------Lasso---------#
	lasso.cv <- cv.glmnet(as.matrix(scz.bio.use),abs(target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min = as.matrix(coef(lasso.cv, s = "lambda.min"))
	coef_list.min <- unlist(coef_list.min)[-1]
	z.combine.lasso.min.bio=scz.bio.use[,coef_list.min!=0]
#assessment-fdr.reg-bio
	 features.se=SEfromHessian(fdr.target.theore.bio$model$hessian)
	 features.coef=	fdr.target.theore.bio$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment2 = data.frame(colnames(z.combine.bio),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment2) = c('feature','p2','beta2','se2')
	 assessment_sum = merge(assessment,assessment2,by='feature',all=T)

#completion of analysis
	target.gene.order$FDR.the = fdr.target.theore$FDR
	target.gene.order$FDR.emp = fdr.target.empiri$FDR
	target.gene.order$qval = p.adjust(target.gene.order$P,method="fdr",n=length(target.gene.order$P))
	target.gene.order$bio.FDR.the = fdr.target.theore.bio$FDR
	target.gene.order$bio.FDR.emp = fdr.target.empiri.bio$FDR

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
	if (!is.null(intersect(colnames(scz.bio.use),colnames(z.combine.lasso.min.bio)))){
	z.combine.bio.lasso=as.matrix(cbind(abs(z.combine),abs(z.combine.lasso.min.bio)))
	set.seed(25)
	fdr.target.theore.bio.lasso=FDRreg(target.z,z.combine.bio.lasso,nulltype = 'theoretical',method = 'pr')
	fdr.target.empiri.bio.lasso=FDRreg(target.z,z.combine.bio.lasso,nulltype = 'empirical',method = 'pr')
	#assessment-fdr.reg-lasso
	 features.se=SEfromHessian(fdr.target.theore.bio.lasso$model$hessian)
	 features.coef=fdr.target.theore.bio.lasso$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment3 = data.frame(colnames(z.combine.bio.lasso),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 colnames(assessment3) = c('feature','p3','beta3','se3')
	 assessment_sum = merge(assessment_sum,assessment3,by='feature',all=T)
	 target.gene.order$bio.FDR.the.lasso = fdr.target.theore.bio.lasso$FDR
	 target.gene.order$bio.FDR.emp.lasso = fdr.target.empiri.bio.lasso$FDR
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

#Save result	
	fwrite(as.data.frame(target.gene.order),'bd2018.gene.bio.fdrreg.txt',sep=' ')
	fwrite(as.data.frame(assessment_sum),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/magma-bd2018.csv',sep=',')
	sum_data = data.frame(qval,fdr.the,fdr.emp,bio.fdr.the,bio.fdr.emp,bio.fdr.the.lasso,bio.fdr.emp.lasso)
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/magma_bd2018.csv',sep=',')	
	save(fdr.target.theore,fdr.target.empiri,fdr.target.theore.bio,fdr.target.empiri.bio,fdr.target.theore.bio.lasso,fdr.target.empiri.bio.lasso,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/magma_bd2018.RData")
#exit the software
	q('no') 