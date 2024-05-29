#####---------Run the FDR-regression---------#####
#R-353
	library(FDRreg);library(data.table);library(powerplus);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/mdd2013')

#----Loading data----#
#Target
	target=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2013/mdd2013.overlap.txt')

#Variants with no sample overlap
	listnames_no = c('anx','asd2019','cannabis','dx.dep','ever.csh','ever.sh','ever.st','insom2019','inte','life.threat.accident','neuroti2018','no.dep','physical.crime','psycho.experi','ptsd2019','scz.eas2019','st','st.dep','wal.scz2018')
	files_no=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2013/',listnames_no,'.overlap.txt'),fread,stringsAsFactor=F)
#Variants with sample overlap
	listnames_with = c('adhd2019','alco2018','antisocial','bd2018','chrono','ed2019','neo.c')
	files_with =lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2013/',listnames_with,'.overlap.txt'),fread,stringsAsFactor=F)
	

#----Get intercept from different comparisons in LDSC program(de-correlation and target only)----#
	mdd2013.inter=c(1,0.029,8.78E-02,0.0128,1.26E-01,-0.0147,-0.0399,-0.0914)
	adhd2019.inter=c(0.029,1,0.0267,-0.0001,0.0337,0.0048,-0.1103,0.0062)
	alco2018.inter=c(0.08777,0.0267,1,0.019414,0.05083,0.0083471,-0.021109,-0.030821)
	antisocial.inter=c(0.0128,-0.0001,1.94E-02,1,-0.0022,0.0061,-0.0019,-0.0061)
	bd2018.inter=c(0.12642,0.0337,5.08E-02,-0.0022,1,-0.00080075,-0.049531,-0.0027653)
	chrono.inter=c(-0.0147,0.0048,8.35E-03,0.0061,-8.01E-04,1,0.0126,0.0012)
	ed2019.inter=c(-0.0399,-0.1103,-2.11E-02,-0.0019,-4.95E-02,0.0126,1,0.0059)
	neo.c.inter=c(-0.0914,0.0062,-3.08E-02,-0.0061,-2.77E-03,0.0012,0.0059,1)


#----Construction of cor.matrix based on those intercept above----#
	cor.matr=matrix(c(mdd2013.inter,adhd2019.inter,alco2018.inter,antisocial.inter,bd2018.inter,chrono.inter,ed2019.inter,neo.c.inter),nrow=8,byrow=T)

#----Calculate de-correlated z.stat based on Matpow function----#
	z.stat.matr=cbind(target$z,files_with[[1]]$z)
	for (i in 2:6){
	z.stat.matr=cbind(z.stat.matr,files_with[[i]]$z)
	}
	z.stat.matr = matrix(cbind(z.stat.matr,files_with[[7]]$z),nrow=8,byrow=T)
	z.decor=Matpow(cor.matr,-0.5) %*% z.stat.matr
	
#----Run fdr.reg function(For all data)----#
	#Remember to add the non-overlap data
	features=cbind(abs(t(z.decor[2:8,])))
	for (i in 1:length(files_no)){features=cbind(features,files_no[[i]]$z)}
	z.decor.target.z=cbind(z.decor[1,])

	colnames(features) = c('adhd2019','alco2018','antisocial','bd2018','chrono','ed2019','neo.c','anx','asd2019','cannabis','dx.dep','ever.csh','ever.sh','ever.st','insom2019','inte','life.threat.accident','neuroti2018','no.dep','physical.crime','psycho.experi','ptsd2019','scz.eas2019','st','st.dep','wal.scz2018')

# #-----lasso for chosing the best variables----#
	lasso.cv.bio <- cv.glmnet(features,abs(z.decor.target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
	coef_list.min.bio <- unlist(coef_list.min.bio)[-1] 
	z.combine.lasso=features[,coef_list.min.bio != 0]
	
	cat('Variables not selected from gwas:',setdiff(colnames(features),colnames(z.combine.lasso)),'\n')
	set.seed(14)
	fdr.st=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'theoretical',method = 'pr')
	fdr.st2=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'empirical',method = 'pr')
	#assessment-fdrr
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.st$model$hessian)
	 features.coef=fdr.st$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 cat(features.pvalue,seq='\n')
	 assessment = cbind(colnames(features),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 fwrite(as.data.frame(assessment),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/snp-mdd2013.csv',sep=',')
#----Get the qval----#
	target$qval = p.adjust(target$pval,method="fdr",n=length(target$pval))
	target$FDR.the = fdr.st$FDR
	target$FDR.emp = fdr.st2$FDR
	fwrite(as.data.frame(target),'mdd2013.fdrreg.txt',sep=' ')

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
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/snp_mdd2013.csv',sep=',')
	
	save(fdr.st,fdr.st2,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/snp_mdd2013.RData")

#----Get the magma data----#
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013')
	
	#Target
	target.decor.z=cbind(z.decor[1,],target)
	colnames(target.decor.z)[1]='z.decor'
	target.decor.z$p.decor=2*pnorm(abs(target.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(target.decor.z),'mdd2013.overlap.4magma.txt',sep=' ')	
	
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