#####---------Run the FDR-regression---------#####
#R-353
	library(FDRreg);library(data.table);library(powerplus);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2012')

#----Loading data----#
#Target
	target=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz2012/scz2012.overlap.txt')

#Variants with no sample overlap
	listnames_no = c('asd2019','cannabis','chrono','dx.dep','edu','ever.csh','ever.sh','ever.st','inte','life.threat.accident','mddco','neo.o','neuroti2018','no.dep','physical.crime','ptsd2019','sleep.dura','st.dep')
	files_no=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz2012/',listnames_no,'.overlap.txt'),fread,stringsAsFactor=F)
#Variants with sample overlap
	listnames_with = c('adhd2019','alco2018','bd2018','ed2019','mdd2019','ocd2018','psycho.experi','reac','tense')
	files_with =lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz2012/',listnames_with,'.overlap.txt'),fread,stringsAsFactor=F)
	

#----Get intercept from different comparisons in LDSC program(de-correlation and target only)----#
	scz2012.inter=c(1,0.0175,0.036717,0.17252,-0.0392,0.0374,0.0323,-0.0002,-0.024,-0.0224)
	adhd2019.inter=c(0.0175,1,0.0267,0.0337,-0.1103,0.0756,0.0075,0.0071,-0.0014,-0.0131)
	alco2018.inter=c(0.036717,0.0267,1,0.05083,-0.021109,0.013106,0.025988,-0.0041905,-0.01252,0.0010582)
	bd2018.inter=c(0.17252,0.0337,0.05083,1,-0.049531,0.050534,0.038071,0.0065843,-0.0065881,-0.027756)
	ed2019.inter=c(-0.0392,-0.1103,-0.021109,-0.049531,1,-0.0621,-0.0109,-0.0089,-0.0091,0.0292)
	mdd2019.inter=c(0.0374,0.0756,0.013106,0.050534,-0.0621,1,0.0087,0.0407,-0.0069,-0.1424)
	ocd2018.inter=c(0.0323,0.0075,0.025988,0.038071,-0.0109,0.0087,1,-0.0025,0.0111,-0.0092)
	psycho.experi.inter=c(-0.0219,0.0071,-0.0041905,0.0065843,-0.0089,0.0407,-0.0025,1,0.0069,-0.0525)
	reac.inter=c(-0.0219,-0.0014,-0.01252,-0.0065881,-0.0091,-0.0069,0.0111,0.0069,1,-0.0525)
	tense.inter=c(-0.0224,-0.0131,0.0010582,-0.027756,0.0292,-0.1424,-0.0092,-0.0167,0.0069,1)

#----Construction of cor.matrix based on those intercept above----#
	cor.matr=matrix(c(scz2012.inter,adhd2019.inter,alco2018.inter,bd2018.inter,ed2019.inter,mdd2019.inter,ocd2018.inter,psycho.experi.inter,reac.inter,tense.inter),nrow=10,byrow=T)

#----Calculate de-correlated z.stat based on Matpow function----#
	z.stat.matr=cbind(target$z,files_with[[1]]$z)
	for (i in 2:8){
	z.stat.matr=cbind(z.stat.matr,files_with[[i]]$z)
	}
	z.stat.matr = matrix(cbind(z.stat.matr,files_with[[9]]$z),nrow=10,byrow=T)
	z.decor=Matpow(cor.matr,-0.5) %*% z.stat.matr
	
#----Run fdr.reg function(For all data)----#
	#Remember to add the non-overlap data
	features=cbind(abs(t(z.decor[2:10,])))
	for (i in 1:length(files_no)){features=cbind(features,files_no[[i]]$z)}
	z.decor.target.z=cbind(z.decor[1,])

	colnames(features) = c('asd2019','cannabis','chrono','dx.dep','edu','ever.csh','ever.sh','ever.st','inte','life.threat.accident','mddco','neo.o','neuroti2018','no.dep','physical.crime','ptsd2019','sleep.dura','st.dep','adhd2019','alco2018','bd2018','ed2019','mdd2019','ocd2018','psycho.experi','reac','tense')

#-----lasso for chosing the best variables----#
	lasso.cv.bio <- cv.glmnet(features,abs(z.decor.target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
	coef_list.min.bio <- unlist(coef_list.min.bio)[-1] 
	z.combine.lasso=features[,coef_list.min.bio != 0]
	
	cat('Variables not selected from gwas:',setdiff(colnames(features),colnames(z.combine.lasso)),'\n')
	set.seed(16)
	fdr.st=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'theoretical',method = 'pr')
	fdr.st2=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'empirical',method = 'pr')
	#assessment-fdrr
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.st$model$hessian)
	 features.coef=fdr.st$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 assessment = cbind(colnames(features),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 fwrite(as.data.frame(assessment),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/snp-scz2012.csv',sep=',')	

#----Get the qval----#
	target$qval = p.adjust(target$pval,method="fdr",n=length(target$pval))
	target$FDR.the = fdr.st$FDR
	target$FDR.emp = fdr.st2$FDR
	fwrite(as.data.frame(target),'scz2012.fdrreg.txt',sep=' ')

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
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/snp_scz2012.csv',sep=',')
	
	save(fdr.st,fdr.st2,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/snp_scz2012.RData")

#----Get the magma data----#
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2012')
	
	#Target
	target.decor.z=cbind(z.decor[1,],target)
	colnames(target.decor.z)[1]='z.decor'
	target.decor.z$p.decor=2*pnorm(abs(target.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(target.decor.z),'scz2012.overlap.4magma.txt',sep=' ')	
	
	#Variants with sample overlap
	for (i in 1:length(files_with)){
	with.decor.z=cbind(z.decor[i+1,],files_with[[i]])
	colnames(with.decor.z)[1]='z.decor'
	name = listnames_with[i]
	with.decor.z$p.decor=2*pnorm(abs(with.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(with.decor.z),paste0(name,'.overlap.4magma.txt'),sep=' ')
	}
	
	#Variants with no sample overlap
	for (i in 1:length(files_no)){
	no.decor.z=files_no[[i]]
	name = listnames_no[i]
	no.decor.z$p.decor=no.decor.z$pval
	no.decor.z$z.decor=no.decor.z$z
    fwrite(as.data.frame(no.decor.z),paste0(name,'.overlap.4magma.txt'),sep=' ')
	}
	
	q('no') 