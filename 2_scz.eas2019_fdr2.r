#####---------Run the FDR-regression---------#####
#R-353
	library(FDRreg);library(data.table);library(powerplus);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz.eas2019')

#----Loading data----#
#Target
	target=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz.eas2019/scz.eas2019.overlap.txt')

#Variants with no sample overlap
	listnames_no = c('alco2018','asd2019','bd2018','cannabis','chrono','ed2019','edu','insom2019','inte','life.threat.accident','mdd2019','mddco','neuroti2018','no.dep','ocd2018','physical.crime','ptsd2019','reac','sleep.dura')
	files_no=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz.eas2019/',listnames_no,'.overlap.txt'),fread,stringsAsFactor=F)
#Variants with sample overlap
	listnames_with = c('ever.st')
	files_with =lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz.eas2019/',listnames_with,'.overlap.txt'),fread,stringsAsFactor=F)
	

#----Get intercept from different comparisons in LDSC program(de-correlation and target only)----#
	scz.eas2019.inter=c(1,0.010536928)
	ever.st.inter=c(0.010536928,1)

#----Construction of cor.matrix based on those intercept above----#
	cor.matr=matrix(c(scz.eas2019.inter,ever.st.inter),nrow=2,byrow=T)

#construction of z.stat.matrix based on gwas summary results#
	#for (i in c(1:12,14:18)){files_no[[i]]$z.stat=files_no[[i]]$beta/files_no[[i]]$se}
	#for (i in c(1,3:8)){files_with[[i]]$z.stat=files_with[[i]]$beta/files_with[[i]]$se}
	#files_no[[13]]$z.stat = files_no[[13]]$z
	#files_with[[2]]$z.stat = files_with[[2]]$z
	#target$z.stat = target$beta/target$se

#----Calculate de-correlated z.stat based on Matpow function----#
	z.stat.matr=matrix(cbind(target$z,files_with[[1]]$z),nrow=2,byrow=T)
	z.decor=Matpow(cor.matr,-0.5) %*% z.stat.matr
	
#----Run fdr.reg function(For all data)----#
	#Remember to add the non-overlap data
	features=abs(z.decor[2,])
	for (i in 1:length(files_no)){features=cbind(features,files_no[[i]]$z)}
	z.decor.target.z=cbind(z.decor[1,])

	colnames(features) = c('ever.st','alco2018','asd2019','bd2018','cannabis','chrono','ed2019','edu','insom2019','inte','life.threat.accident','mdd2019','mddco','neuroti2018','no.dep','ocd2018','physical.crime','ptsd2019','reac','sleep.dura')

#-----lasso for chosing the best variables----#
	lasso.cv.bio <- cv.glmnet(features,abs(z.decor.target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
	coef_list.min.bio <- unlist(coef_list.min.bio)[-1] 
	z.combine.lasso=features[,coef_list.min.bio != 0]
	
	cat('Variables not selected from gwas:',setdiff(colnames(features),colnames(z.combine.lasso)),'\n')
	set.seed(18)
	fdr.st=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'theoretical',method = 'pr')
	fdr.st2=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'empirical',method = 'pr')
	#assessment-fdrr
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.st$model$hessian)
	 features.coef=fdr.st$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment = cbind(colnames(features),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 fwrite(as.data.frame(assessment),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/snp-scz.eas2019.csv',sep=',')	

#----Get the qval----#
	target$qval = p.adjust(target$pval,method="fdr",n=length(target$pval))
	target$FDR.the = fdr.st$FDR
	target$FDR.emp = fdr.st2$FDR
	fwrite(as.data.frame(target),'scz.eas2019.fdrreg.txt',sep=' ')

#----Compare FDR-reg result----#
	qval <- c(length(which(target$qval<0.5)),
	length(which(target$qval<0.4)),
	length(which(target$qval<0.3)),
	length(which(target$qval<0.2)),
	length(which(target$qval<0.1)),
	length(which(target$qval<0.05)),
	length(which(target$qval<0.04)),
	length(which(target$qval<0.03)),
	length(which(target$qval<0.02)),
	length(which(target$qval<0.01)),
	length(which(target$qval<0.001)),
	length(which(target$qval<5e-04)),
	length(which(target$qval<5e-06)),
	length(which(target$qval<5e-08)))
	fdr.the <- c(length(which(target$FDR.the<0.5)),
	length(which(target$FDR.the<0.4)),
	length(which(target$FDR.the<0.3)),
	length(which(target$FDR.the<0.2)),
	length(which(target$FDR.the<0.1)),
	length(which(target$FDR.the<0.05)),
	length(which(target$FDR.the<0.04)),
	length(which(target$FDR.the<0.03)),
	length(which(target$FDR.the<0.02)),
	length(which(target$FDR.the<0.01)),
	length(which(target$FDR.the<0.001)),
	length(which(target$FDR.the<5e-04)),
	length(which(target$FDR.the<5e-06)),
	length(which(target$FDR.the<5e-08)))
	fdr.emp <- c(length(which(target$FDR.emp<0.5)),
	length(which(target$FDR.emp<0.4)),
	length(which(target$FDR.emp<0.3)),
	length(which(target$FDR.emp<0.2)),
	length(which(target$FDR.emp<0.1)),
	length(which(target$FDR.emp<0.05)),
	length(which(target$FDR.emp<0.04)),
	length(which(target$FDR.emp<0.03)),
	length(which(target$FDR.emp<0.02)),
	length(which(target$FDR.emp<0.01)),
	length(which(target$FDR.emp<0.001)),
	length(which(target$FDR.emp<5e-04)),
	length(which(target$FDR.emp<5e-06)),
	length(which(target$FDR.emp<5e-08)))
	sum_data = data.frame(qval,fdr.the,fdr.emp)
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/snp_scz.eas2019.csv',sep=',')
	
	save(fdr.st,fdr.st2,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/snp_scz.eas2019.RData")

#----Get the magma data----#
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz.eas2019')
	
	#Target
	target.decor.z=cbind(z.decor[1,],target)
	colnames(target.decor.z)[1]='z.decor'
	target.decor.z$p.decor=2*pnorm(abs(target.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(target.decor.z),'scz.eas2019.overlap.4magma.txt',sep=' ')	
	
	#Variants with sample overlap
	with.decor.z=cbind(z.decor[2,],files_with[[1]])
	colnames(with.decor.z)[1]='z.decor'
	name = listnames_with[1]
	with.decor.z$p.decor=2*pnorm(abs(with.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(with.decor.z),paste0(name,'.overlap.4magma.txt'),sep=' ')
	
	
	#Variants with no sample overlap
	for (i in 1:length(files_no)){
	no.decor.z=files_no[[i]]
	name = listnames_no[i]
	no.decor.z$p.decor=no.decor.z$pval
	no.decor.z$z.decor=no.decor.z$z
    fwrite(as.data.frame(no.decor.z),paste0(name,'.overlap.4magma.txt'),sep=' ')
	}
	
	q('no') 