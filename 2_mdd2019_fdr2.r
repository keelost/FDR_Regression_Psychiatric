#####---------Run the FDR-regression---------#####
#R-353
	library(FDRreg);library(data.table);library(powerplus);library(glmnet)
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/mdd2019')

#----Loading data----#
#Target
	target=fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2019/mdd2019.overlap.txt')

#Variants with no sample overlap
	listnames_no = c('antisocial','college2013','edu','inte','life.threat.accident','ocd2018','sa.ipsych','scz.eas2019','wit.death')
	files_no=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2019/',listnames_no,'.overlap.txt'),fread,stringsAsFactor=F)
#Variants with sample overlap
	listnames_with = c('ad','adhd2019','alco2018','anx','asd2019','bd2018','cannabis','chrono','dx.dep','dx.life.threat.ill','ed2019','ever.csh','ever.sh','ever.st','insom2019','neuroti2018','no.dep','physical.crime','psycho.experi','ptsd2019','sleep.dura','st','st.dep','tense','wal.scz2018')
	files_with =lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2019/',listnames_with,'.overlap.txt'),fread,stringsAsFactor=F)

#----Get intercept from different comparisons in LDSC program(de-correlation and target only)----#
	mdd2019.inter=c(1,0.0174,0.0756,1.31E-02,-0.1296,0.0511,5.05E-02,0.0179,0.0393,0.2599,1.49E-02,-0.0621,0.122,0.0811,0.1404,0.085598,0.2712,0.0454,0.0282,0.0407,0.0488,-0.0148,0.064,0.0295,-0.1424,0.0361)
	ad.inter=c(0.0174,1,0.0022,-0.0043,-0.0118,0.0029,-0.0019,-0.005,0.0029,0.0071,0.0104,-0.0038,0.0059,0.0097,0.0131,0.0078,0.023,-0.0009,0.0066,0.0039,0.0121,-0.0078,-0.0004,-0.0627,-0.0525,-0.0423)
	adhd2019.inter=c(0.0756,0.0022,1,0.0267,-0.0015,0.3392,0.0337,-0.0106,0.0048,0.012,-0.0033,-0.1103,0.0089,0.0165,0.0046,0.0049,0.0093,-0.0057,0.002,0.0071,0.0258,-0.0032,0.0062,0.0075,-0.0131,0.0226)
	alco2018.inter=c(0.013106,-0.0043,0.0267,1,-0.0031464,0.0016723,0.05083,0.024314,0.0083471,-0.0024438,0.00066851,-0.021109,0.0046294,0.0023061,0.0074556,0.0029023,0.0031389,-0.00047954,-0.0047122,-0.0041905,0.019928,-0.0006848,0.0031995,0.0035842,0.0010582,0.037279)
	anx.inter=c(-0.1296,-0.0118,-0.0015,-3.15E-03,1,0.0001,-0.0147,0.0076,-0.008,-0.0651,-0.0067,0.0255,-0.0412,-0.0299,-0.057,-0.1019,-0.4392,-0.03,-0.0107,-0.0133,-0.0275,0.0162,-0.034,-0.022,0.3245,-0.0141)
	asd2019.inter=c(0.0511,0.0029,0.3392,1.67E-03,0.0001,1,0.0103,-0.0015,-0.0019,0.003,0.0104,-0.0925,0.0038,0.0051,-0.0023,0.0007,-0.0018,-0.0028,0.0047,0.0001,0.023,0.0008,-0.0003,0.0072,-0.0077,0.0112)
	bd2018.inter=c(0.050534,-0.0019,0.0337,5.08E-02,-0.0147,0.0103,1,0.0025462,-0.00080075,0.011231,0.0064019,-0.049531,0.0033674,0.000097317,0.010365,-0.01329,0.010506,-0.0048065,0.0054808,0.0065843,0.022154,-0.004115,0.0093393,0.0043984,-0.027756,0.17745)
	cannabis.inter=c(0.0179,-0.005,-0.0106,2.43E-02,0.0076,-0.0015,2.55E-03,1,0.0444,0.0445,0.0221,0.0205,0.0934,0.0739,0.0929,-0.0117,-0.0108,0.0153,0.103,0.0367,0.0017,0.0067,0.0239,0.0092,0.0012,0.0066)
	chrono.inter=c(0.0393,0.0029,0.0048,8.35E-03,-0.008,-0.0019,-8.01E-04,0.0444,1,0.0273,0.0104,0.0126,0.0246,0.0153,0.0264,0.0081,0.0739,0.03,0.0274,0.0048,0.0055,0.0197,0.0219,0.0164,-0.0236,-0.004)
	dx.dep.inter=c(0.2599,0.0071,0.012,-2.44E-03,-0.0651,0.003,1.12E-02,0.0445,0.0273,1,0.0301,-0.0254,0.3335,0.21,0.3698,0.0438,0.1522,0.056,0.0393,0.0831,0.0716,-0.0008,0.1794,0.0877,-0.065,0.0079)
	dx.life.threat.ill.inter=c(0.014924,0.0104,-0.0033,6.69E-04,-0.0067,0.0104,6.40E-03,0.0221,0.0104,0.0301,1,-0.0071401,0.028958,0.021177,0.042984,0.014742,0.010186,-0.0028057,0.034284,0.032196,0.021377,-0.0035839,0.018483,0.024483,-0.007445,0.00015552)
	ed2019.inter=c(-0.0621,-0.0038,-0.1103,-2.11E-02,0.0255,-0.0925,-4.95E-02,0.0205,0.0126,-0.0254,-7.14E-03,1,-0.024,-0.0123,-0.0177,-0.0014,-0.0173,-0.006,-0.0042,-0.0089,-0.0357,-0.0016,-0.0233,-0.0071,0.0292,-0.03)
	ever.csh.inter=c(0.122,0.0059,0.0089,4.63E-03,-0.0412,0.0038,3.37E-03,0.0934,0.0246,0.3335,2.90E-02,-0.024,1,0.4957,0.5879,0.043,0.118,0.059,0.0729,0.1221,0.0772,-0.0026,0.3438,0.2613,-0.0515,0.009)
	ever.sh.inter=c(0.0811,0.0097,0.0165,2.31E-03,-0.0299,0.0051,9.73E-05,0.0739,0.0153,0.21,2.12E-02,-0.0123,0.4957,1,0.285,0.0253,0.0672,0.0352,0.051,0.0802,0.0455,0.0007,0.1768,0.1312,-0.0307,0.0088)
	ever.st.inter=c(0.1404,0.0131,0.0046,7.46E-03,-0.057,-0.0023,1.04E-02,0.0929,0.0264,0.3698,4.30E-02,-0.0177,0.5879,0.285,1,0.0516,0.1637,0.0634,0.0764,0.1195,0.0828,-0.005,0.3506,0.2812,-0.0605,0.0154)
	insom2019.inter=c(0.085598,0.0078,0.0049,2.90E-03,-0.1019,0.0007,-1.33E-02,-0.0117,0.0081,0.0438,1.47E-02,-0.0014,0.043,0.0253,0.0516,1,0.21271,0.02679,0.0046977,0.014895,0.017123,-0.116,0.051445,0.024601,-0.089405,0.00043075)
	neuroti2018.inter=c(0.2712,0.023,0.0093,3.14E-03,-0.4392,-0.0018,1.05E-02,-0.0108,0.0739,0.1522,1.02E-02,-0.0173,0.118,0.0672,0.1637,0.21271,1,0.0893,0.027,0.0322,0.041,-0.0259,0.0941,0.0596,-0.4012,0.005)
	no.dep.inter=c(0.0454,-0.0009,-0.0057,-4.80E-04,-0.03,-0.0028,-4.81E-03,0.0153,0.03,0.056,-2.81E-03,-0.006,0.059,0.0352,0.0634,0.02679,0.0893,1,0.016,0.0148,0.0076,0.0005,0.0282,0.0312,-0.0414,-0.0029)
	physical.crime.inter=c(0.0282,0.0066,0.002,-4.71E-03,-0.0107,0.0047,5.48E-03,0.103,0.0274,0.0393,3.43E-02,-0.0042,0.0729,0.051,0.0764,0.0046977,0.027,0.016,1,-8.0627,0.0188,-0.004,0.0393,0.021,-0.0228,0.0128)
	psycho.experi.inter=c(0.0407,0.0039,0.0071,-4.19E-03,-0.0133,0.0001,6.58E-03,0.0367,0.0048,0.0831,3.22E-02,-0.0089,0.1221,0.0802,0.1195,0.014895,0.0322,0.0148,0.0468,1,0.0398,-0.0006,0.0931,-0.0627,-0.0525,-0.0423)
	ptsd2019.inter=c(0.0488,0.0121,0.0258,1.99E-02,-0.0275,0.023,2.22E-02,0.0017,0.0055,0.0716,2.14E-02,-0.0357,0.0772,0.0455,0.0828,0.017123,0.041,0.0076,0.0188,0.0398,1,-0.0028,0.0854,0.0468,-0.0294,0.0198)
	sleep.dura.inter=c(-0.0148,-0.0078,-0.0032,-6.85E-04,0.0162,0.0008,-4.12E-03,0.0067,0.0197,-0.0008,-3.58E-03,-0.0016,-0.0026,0.0007,-0.005,-0.116,-0.0259,0.0005,-0.004,-0.0006,-0.0028,1,-0.0059,0.002,0.0017,-0.0135)
	st.inter=c(0.064,-0.0004,0.0062,3.20E-03,-0.034,-0.0003,9.34E-03,0.0239,0.0219,0.1794,1.85E-02,-0.0233,0.3438,0.1768,0.3506,0.051445,0.0941,0.0282,0.0393,0.0931,0.0854,-0.0059,1,0.1765,-0.0361,-0.0423)
	st.dep.inter=c(0.0295,-0.0014,0.0075,3.58E-03,-0.022,0.0072,4.40E-03,0.0092,0.0164,0.0877,2.45E-02,-0.0071,0.2613,0.1312,0.2812,0.024601,0.0596,0.0312,0.021,0.0567,0.0468,0.002,0.1765,1,-0.021424,-0.0423)
	tense.inter=c(-0.1424,-0.0186,-0.0131,1.06E-03,0.3245,-0.0077,-2.78E-02,0.0012,-0.0236,-0.065,-7.45E-03,0.0292,-0.0515,-0.0307,-0.0605,-0.089405,-0.4012,-0.0414,-0.0228,-0.0167,-0.0294,0.0017,-0.0361,-0.021424,1,-0.0423)
	wal.scz2018.inter=c(0.0361,0.0123,0.0226,3.73E-02,-0.0141,0.0112,1.77E-01,0.0066,-0.004,0.0079,1.56E-04,-0.03,0.009,0.0088,0.0154,0.00043075,0.005,-0.0029,0.0128,0.0069,0.0198,-0.0135,0.0044,-0.0011587,-0.041,1)


#----Construction of cor.matrix based on those intercept above----#
	cor.matr=matrix(c(mdd2019.inter,ad.inter,adhd2019.inter,alco2018.inter,anx.inter,asd2019.inter,bd2018.inter,cannabis.inter,chrono.inter,dx.dep.inter,dx.life.threat.ill.inter,ed2019.inter,ever.csh.inter,ever.sh.inter,ever.st.inter,insom2019.inter,neuroti2018.inter,no.dep.inter,physical.crime.inter,psycho.experi.inter,ptsd2019.inter,sleep.dura.inter,st.inter,st.dep.inter,tense.inter,wal.scz2018.inter),nrow=26,byrow=T)

#----Calculate de-correlated z.stat based on Matpow function----#
	z.stat.matr=cbind(target$z,files_with[[1]]$z)
	for (i in 2:24){
	z.stat.matr=cbind(z.stat.matr,files_with[[i]]$z)
	}
	z.stat.matr = matrix(cbind(z.stat.matr,files_with[[25]]$z),nrow=26,byrow=T)
	z.decor=Matpow(cor.matr,-0.5) %*% z.stat.matr
	
#----Run fdr.reg function(For all data)----#
	#Remember to add the non-overlap data
	features=cbind(abs(t(z.decor[2:26,])))
	for (i in 1:length(files_no)){features=cbind(features,files_no[[i]]$z)}
	z.decor.target.z=cbind(z.decor[1,])

	colnames(features) = c('ad','adhd2019','alco2018','anx','asd2019','bd2018','cannabis','chrono','dx.dep','dx.life.threat.ill','ed2019','ever.csh','ever.sh','ever.st','insom2019','neuroti2018','no.dep','physical.crime','psycho.experi','ptsd2019','sleep.dura','st','st.dep','tense','wal.scz2018','antisocial','college2013','edu','inte','life.threat.accident','ocd2018','sa.ipsych','scz.eas2019','wit.death')

#-----lasso for chosing the best variables----#
	lasso.cv.bio <- cv.glmnet(features,abs(z.decor.target.z),family='gaussian',nlambda=50,alpha=1,standardize=T,parallel=T)
	coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
	coef_list.min.bio <- unlist(coef_list.min.bio)[-1] 
	z.combine.lasso=features[,coef_list.min.bio != 0]
	
	cat('Variables not selected from gwas:',setdiff(colnames(features),colnames(z.combine.lasso)),'\n')
	set.seed(15)
	fdr.st=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'theoretical',method = 'pr')
	fdr.st2=FDRreg(z.decor.target.z,z.combine.lasso,nulltype = 'empirical',method = 'pr')
	#assessment-fdrr
	 library(HelpersMG)
	 features.se=SEfromHessian(fdr.st$model$hessian)
	 features.coef=fdr.st$model$coef
	 features.z.score=features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	 features.pvalue=2*pnorm(abs(features.z.score),lower.tail=F)	  
	 assessment = cbind(colnames(features),features.pvalue,features.coef[c(2:length(features.coef))],features.se[c(2:length(features.se))])
	 fwrite(as.data.frame(assessment),'/exeh_4/jinghong_qiu/SO_Lab/Contribution/snp-mdd2019.csv',sep=',')
#----Get the qval----#
	target$qval = p.adjust(target$pval,method="fdr",n=length(target$pval))
	target$FDR.the = fdr.st$FDR
	target$FDR.emp = fdr.st2$FDR
	fwrite(as.data.frame(target),'mdd2019.fdrreg.txt',sep=' ')

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
	fwrite(as.data.frame(sum_data),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/snp_mdd2019.csv',sep=',')
	
	save(fdr.st,fdr.st2,file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/snp_mdd2019.RData")

#----Get the magma data----#
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2019')
	
	#Target
	target.decor.z=cbind(z.decor[1,],target)
	colnames(target.decor.z)[1]='z.decor'
	target.decor.z$p.decor=2*pnorm(abs(target.decor.z$z.decor),lower.tail=F)
    fwrite(as.data.frame(target.decor.z),'mdd2019.overlap.4magma.txt',sep=' ')	
	
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