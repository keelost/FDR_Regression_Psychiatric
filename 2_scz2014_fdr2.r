#####---------Run the FDR-regression---------#####
#R-353
	library(FDRreg);library(data.table);library(powerplus);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2014')

#----Loading data----#
#Target
	target=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz2014/scz2014.overlap.txt')

#Variants with no sample overlap
	listnames_no = c('ad','asd2019','cannabis','chrono','dx.dep','dx.life.threat.ill','edu','ever.csh','ever.sh','inte','mddco','neuroti2018','no.dep','pd','sleep.dura','smoking','st','st.dep','wit.death')
	files_no=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz2014/',listnames_no,'.overlap.txt'),fread,stringsAsFactor=F)
#Variants with sample overlap
	listnames_with = c('adhd2019','alco2018','anx','bd2018','ed2019','ever.st','life.threat.accident','mdd2019','ocd2018','physical.crime','psycho.experi','ptsd2019','reac','tense')
	files_with =lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/scz2014/',listnames_with,'.overlap.txt'),fread,stringsAsFactor=F)
	

#----Get intercept from different comparisons in LDSC program(de-correlation and target only)----#
	scz2014.inter=c(1,0.0349,0.036723,-0.0147,0.21901,-0.0348,0.0127,0.011,0.0424,0.0294,0.0163,0.0049,0.0216,-0.0156,-0.0404)
	adhd2019.inter=c(0.0349,1,0.0267,-0.0015,0.0337,-0.1103,0.0046,0.0075,0.0756,0.0075,0.002,0.0071,0.0258,-0.0014,-0.0131)
	alco2018.inter=c(0.036723,0.0267,1,-0.0031464,0.05083,-0.021109,0.0074556,0.0063964,0.013106,0.025988,-0.0047122,-0.0041905,0.019928,-0.01252,0.0010582)
	anx.inter=c(-0.0147,-0.0015,-0.0031464,1,-0.0147,0.0255,-0.057,0.007,-0.1296,-0.0071,-0.0107,-0.0133,-0.0275,0.0029,0.3245)
	bd2018.inter=c(0.21901,0.0337,0.05083,-0.0147,1,-0.049531,0.010365,0.007435,0.050534,0.038071,0.0054808,0.0065843,0.022154,-0.0065881,-0.027756)
	ed2019.inter=c(-0.0348,-0.1103,-0.021109,0.0255,-0.049531,1,-0.0177,-0.0079,-0.0621,-0.0109,-0.0042,-0.0089,-0.0357,-0.0091,0.0292)
	ever.st.inter=c(0.0127,0.0046,0.0074556,-0.057,0.010365,-0.0177,1,0.0668,0.1404,-0.0026,0.0764,0.1195,0.0828,0.0024,-0.0605)
	life.threat.accident.inter=c(0.011,0.0075,0.0063964,0.007,0.007435,-0.0079,0.0668,1,0.0226,-0.0091,0.1115,0.0468,0.0232,0.006,-0.0211)
	mdd2019.inter=c(0.0424,0.0756,0.013106,-0.1296,0.050534,-0.0621,0.1404,0.0226,1,0.0087,0.0282,0.0407,0.0488,-0.0069,-0.1424)
	ocd2018.inter=c(0.0294,0.0075,0.025988,-0.0071,0.038071,-0.0109,-0.0026,-0.0091,0.0087,1,-0.0055,-0.0025,0.0013,0.0111,-0.0092)
	physical.crime.inter=c(0.0163,0.002,-0.0047122,-0.0107,0.0054808,-0.0042,0.0764,0.1115,0.0282,-0.0055,1,-8.0627,0.0188,-0.0035,-0.0228)
	psycho.experi.inter=c(-0.0117,0.0071,-0.0041905,-0.0133,0.0065843,-0.0089,0.1195,0.0468,0.0407,-0.0025,0.0468,1,0.0398,0.0069,-0.0525)
	ptsd2019.inter=c(0.0216,0.0258,0.019928,-0.0275,0.022154,-0.0357,0.0828,0.0232,0.0488,0.0013,0.0188,0.0398,1,0.0039,-0.0294)
	reac.inter=c(-0.0117,-0.0014,-0.01252,0.0029,-0.0065881,-0.0091,0.0024,0.006,-0.0069,0.0111,-0.0035,0.0069,0.0039,1,-0.0525)
	tense.inter=c(-0.0404,-0.0131,0.0010582,0.3245,-0.027756,0.0292,-0.0605,-0.0211,-0.1424,-0.0092,-0.0228,-0.0167,-0.0294,0.0069,1)

#----Construction of cor.matrix based on those intercept above----#
	cor.matr=matrix(c(scz2014.inter,adhd2019.inter,alco2018.inter,anx.inter,bd2018.inter,ed2019.inter,ever.st.inter,life.threat.accident.inter,mdd2019.inter,ocd2018.inter,physical.crime.inter,psycho.experi.inter,ptsd2019.inter,reac.inter,tense.inter),nrow=15,byrow=T)

#----Calculate de-correlated z.stat based on Matpow function----#
	z.stat.matr=cbind(target$z,files_with[[1]]$z)
	for (i in 2:13){
	z.stat.matr=cbind(z.stat.matr,files_with[[i]]$z)
	}
	z.stat.matr = matrix(cbind(z.stat.matr,files_with[[14]]$z),nrow=15,byrow=T)
	z.decor=Matpow(cor.matr,-0.5) %*% z.stat.matr
	
#----Run fdr.reg function(For all data)----#
	#Remember to add the non-overlap data
	features=cbind(abs(t(z.decor[2:15,])))
	for (i in 1:length(files_no)){features=cbind(features,files_no[[i]]$z)}
	z.decor.target.z=cbind(z.decor[1,])

	colnames(features) = c('adhd2019','alco2018','anx','bd2018','ed2019','ever.st','life.threat.accident','mdd2019','ocd2018','physical.crime','psycho.experi','ptsd2019','reac','tense','ad','asd2019','cannabis','chrono','dx.dep','dx.life.threat.ill','edu','ever.csh','ever.sh','inte','mddco','neuroti2018','no.dep','pd','sleep.dura','smoking','st','st.dep','wit.death')

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
	 fwrite(as.data.frame(assessment),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/snp-scz2014.csv',sep=',')	

#----Get the qval----#
	target$qval = p.adjust(target$pval,method="fdr",n=length(target$pval))
	target$FDR.the = fdr.st$FDR
	target$FDR.emp = fdr.st2$FDR
	fwrite(as.data.frame(target),'scz2014.fdrreg.txt',sep=' ')

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
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/snp_scz2014.csv',sep=',')
	
	save(fdr.st,fdr.st2,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/snp_scz2014.RData")

#----Get the magma data----#
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2014')
	
	#Target
	target.decor.z=cbind(z.decor[1,],target)
	colnames(target.decor.z)[1]='z.decor'
	target.decor.z$p.decor=2*pnorm(abs(target.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(target.decor.z),'scz2014.overlap.4magma.txt',sep=' ')	
	
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