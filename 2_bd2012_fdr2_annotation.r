#####---------Run the FDR-regression---------#####
#R-353
	library(FDRreg);library(data.table);library(powerplus);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/bd2012')

#----Loading data----#
#Target
	target=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2012/bd2012.overlap.txt')

#Variants with no sample overlap
	listnames_no = c('ad','asd2019','cannabis','college2013','dx.dep','edu','ever.csh','ever.sh','ever.st','life.threat.accident','mddco','neuroti2018','no.dep','ocd2018','ptsd2019','scz.eas2019','sleep.dura','wal.scz2018')
	files_bio = fread('/exeh_3/rstao/dr.so/002.meta/002.clear.data/sa.pgc/005.snp.anno/paintor.bd2018.no.1.txt')
	file_bio_fin = files_bio
	for (i in 2:10){
	files_bio = fread(paste0('/exeh_3/rstao/dr.so/002.meta/002.clear.data/sa.pgc/005.snp.anno/paintor.bd2018.no.',i,'.txt'))
	file_bio_fin = rbind(file_bio_fin,files_bio)}
	
	files_no=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2012/',listnames_no,'.overlap.txt'),fread,stringsAsFactor=F)
#Variants with sample overlap
	listnames_with = c('alco2018','ed2019','mdd2019')
	files_with =lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2012/',listnames_with,'.overlap.txt'),fread,stringsAsFactor=F)
	
	
#----Get intercept from different comparisons in LDSC program(de-correlation and target only)----#
	bd2012.inter=c(1,3.07E-02,-0.0429,0.0393)
	alco2018.inter=c(0.030665,1,-0.021109,0.013106)
	ed2019.inter=c(-0.0429,-2.11E-02,1,-0.0621)
	mdd2019.inter=c(0.0393,1.31E-02,-0.0621,1)

#----Construction of cor.matrix based on those intercept above----#
	
	#calculate p-value from LDSC intercept and se, p<0.05 means the two gwas datasets have overlap and are needed to do de-correlation#
	#overlap.z.stat=overlap.intercept/overlap.se #overlap.intercept=LDSC.intercept; overlap.se=LDSC.se in LDSC.intercept#
	#overlap.pval=2*pnorm(abs(overlap.z.stat),lower.tail=F)# 
	cor.matr=matrix(c(bd2012.inter,alco2018.inter,ed2019.inter,mdd2019.inter),nrow=4,byrow=T)
	#head(cor.matr)

#construction of z.stat.matrix based on gwas summary results#
	#for (i in c(1:12,14:18)){files_no[[i]]$z.stat=files_no[[i]]$beta/files_no[[i]]$se}
	#for (i in c(1,3:8)){files_with[[i]]$z.stat=files_with[[i]]$beta/files_with[[i]]$se}
	#files_no[[13]]$z.stat = files_no[[13]]$z
	#files_with[[2]]$z.stat = files_with[[2]]$z
	#target$z.stat = target$beta/target$se

#----Calculate de-correlated z.stat based on Matpow function----#
	z.stat.matr=cbind(target$z,files_with[[1]]$z)
	z.stat.matr=cbind(z.stat.matr,files_with[[2]]$z)
	z.stat.matr = matrix(cbind(z.stat.matr,files_with[[3]]$z),nrow=4,byrow=T)
	z.decor=Matpow(cor.matr,-0.5) %*% z.stat.matr
	
#----Run fdr.reg function(For all data)----#
	#Remember to add the non-overlap data
	features=cbind(abs(t(z.decor[2:4,])))
	for (i in 1:length(files_no)){features=cbind(features,files_no[[i]]$z)}
	z.decor.target.z=cbind(z.decor[1,])

	colnames(features) = c('alco2018','ed2019','mdd2019','ad','asd2019','cannabis','college2013','dx.dep','edu','ever.csh','ever.sh','ever.st','life.threat.accident','mddco','neuroti2018','no.dep','ocd2018','ptsd2019','scz.eas2019','sleep.dura','wal.scz2018')

#-----lasso for chosing the best variables----#
	lasso.cv.bio <- cv.glmnet(file_bio_fin,abs(z.decor.target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
	coef_list.min.bio <- unlist(coef_list.min.bio)[-1] 
	z.combine.lasso=file_bio_fin[,coef_list.min.bio != 0]
	z.combine.final = cbind(features,z.combine.lasso)

	fdr.st=FDRreg(z.decor.target.z,z.combine.final,nulltype = 'theoretical',method = 'pr')
	fdr.st2=FDRreg(z.decor.target.z,z.combine.final,nulltype = 'empirical',method = 'pr')
	
#assessment-fdrr
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.st$model$hessian)
	 features.coef=fdr.st$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)
	 cat(features.pvalue,seq='\n')
	 assessment = cbind(colnames(z.combine.final),features.pvalue)
	 fwrite(as.data.frame(assessment),'/exeh_4/jinghong_qiu/logfile/snp-bd2012.csv',sep=',')
	
#----Get the qval----#
	target$qval = p.adjust(target$pval,method="fdr",n=length(target$pval))
	target$FDR.the = fdr.st$FDR
	target$FDR.emp = fdr.st2$FDR
	fwrite(as.data.frame(target),'bd2012.fdrreg.txt',sep=' ')

#----Compare FDR-reg result----#
	cat(print(paste0('Here is the fdr result of snp\n -----------------------')),sep='\n')
	cat('Here is the number of origin q: ',
	length(which(target$qval<0.5)),
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
	length(which(target$qval<5e-08)),
	'Here is the number of FDR-the: ',
	length(which(target$FDR.the<0.5)),
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
	length(which(target$FDR.the<5e-08)),
	'Here is the number of FDR-emp: ',
	length(which(target$FDR.emp<0.5)),
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
	length(which(target$FDR.emp<5e-08)),
	sep = "\n"
)

# #----Get the magma data----#
	# setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2012')
	
	# #Target
	# target.decor.z=cbind(z.decor[1,],target)
	# colnames(target.decor.z)[1]='z.decor'
	# target.decor.z$p.decor=2*pnorm(abs(target.decor.z$z.decor),lower.tail=F)
    # fwrite(as.data.frame(target.decor.z),'bd2012.overlap.4magma.txt',sep=' ')	
	
	# #Variants with sample overlap
	# for (i in 1:length(files_with)){
	# with.decor.z=cbind(z.decor[i+1,],files_with[[i]])
	# colnames(with.decor.z)[1]='z.decor'
	# name = listnames_with[i]
	# with.decor.z$p.decor=2*pnorm(abs(with.decor.z$z.decor),lower.tail=F)
    # fwrite(as.data.frame(with.decor.z),paste0(name,'.overlap.4magma.txt'),sep=' ')
	# }
	
	# #Variants with no sample overlap
	# for (i in 1:length(files_no)){
	# no.decor.z=files_no[[i]]
	# name = listnames_no[i]
	# no.decor.z$p.decor=no.decor.z$pval
	# no.decor.z$z.decor=no.decor.z$z
    # fwrite(as.data.frame(no.decor.z),paste0(name,'.overlap.4magma.txt'),sep=' ')
	# }
	
	# q('no') 