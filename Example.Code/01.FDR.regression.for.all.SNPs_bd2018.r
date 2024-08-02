#----------------------------------#
#---- FDR-regression for BD2018----#
#----------------------------------#

# SoftWare: R-4.3.3
# Platform: Server 203
# Updated: 2024/08/04
# Author: Jinghong QIU/Shitao RAO/Hon-Cheong SO

#---- load packages and set dir ----#
rm(list = ls())
chooseCRANmirror(ind=25) # for HK mirror setting, getCRANmirrors() can check the number
options(stringsAsFactors = F)
start_time = Sys.time() # time tracking
set.seed(13) # for repeated result
# load packages
library(FDRreg);library(data.table);library(powerplus);library(glmnet)

#---- set dir ----#
setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/bd2018/')

#---- Pipeline ----#
# (1) Load targets
target = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2018/bd2018.overlap.txt')

# Variants with no sample overlap
listnames_no = c('ad', 'asd2019', 'cannabis', 'college2013', 'edu', 'ever.csh', 'ever.sh', 'ever.st', 'life.threat.accident', 'mddco', 'neuroti2018', 'no.dep', 'physical.crime', 'psycho.experi', 'reac', 'scz.eas2019', 'sleep.dura', 'wit.death')
files_no = lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2018/', listnames_no, '.overlap.txt'), fread, stringsAsFactor = F)

# Variants with sample overlap
listnames_with = c('adhd2019','alco2018','dx.dep','ed2019','mdd2019','ocd2018','ptsd2019','tense','wal.scz2018')
files_with = lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2018/', listnames_with, '.overlap.txt'), fread, stringsAsFactor = F)

# (2) Get intercept from different comparisons in LDSC program(de-correlation and target only)
# The intercept is a crucial parameter for detecting and correcting false-positive association signals caused by systematic biases. 
bd2018.inter = c(1,0.0337,5.08E-02,0.011231,-0.049531,0.050534,0.038071,0.022154,-0.027756,0.17745)
adhd2019.inter = c(0.0337,1,0.0267,0.012,-0.1103,0.0756,0.0075,0.0258,-0.0131,0.0226)
alco2018.inter = c(0.05083,0.0267,1,-0.0024438,-0.021109,0.013106,0.025988,0.019928,0.0010582,0.037279)
dx.dep.inter = c(1.12E-02,0.012,-2.44E-03,1,-0.0254,0.2599,0.0031,0.0716,-0.065,0.0079)
ed2019.inter = c(-4.95E-02,-0.1103,-2.11E-02,-0.0254,1,-0.0621,-0.0109,-0.0357,0.0292,-0.03)
mdd2019.inter = c(5.05E-02,0.0756,1.31E-02,0.2599,-0.0621,1,0.0087,0.0488,-0.1424,0.0361)
ocd2018.inter = c(3.81E-02,0.0075,2.60E-02,0.0031,-0.0109,0.0087,1,0.0013,-0.0092,0.0257)
ptsd2019.inter = c(2.22E-02,0.0258,1.99E-02,0.0716,-0.0357,0.0488,0.0013,1,-0.0294,0.0198)
tense.inter = c(-2.78E-02,-0.0131,1.06E-03,-0.065,0.0292,-0.1424,-0.0092,-0.0294,1,-0.0423)
wal.scz2018.inter = c(1.77E-01,0.0226,3.73E-02,0.0079,-0.03,0.0361,0.0257,0.0198,-0.041,1)

# (3) Construction of cor.matrix based on those intercept above
cor.matr = matrix(c(bd2018.inter, adhd2019.inter, alco2018.inter, dx.dep.inter, ed2019.inter, mdd2019.inter, ocd2018.inter, ptsd2019.inter, tense.inter, wal.scz2018.inter), nrow=10, byrow=T)

# (4) Calculate de-correlated z.stat based on Matpow function to avoid the differences caused by sample overlap
z.stat.matr = cbind(target$z, files_with[[1]]$z)

for (i in 2:8){
	z.stat.matr = cbind(z.stat.matr,files_with[[i]]$z)
	}

z.stat.matr = matrix(cbind(z.stat.matr,files_with[[9]]$z), nrow = 10, byrow = T)
z.decor = Matpow(cor.matr,-0.5) %*% z.stat.matr

# (5) Run fdr.reg function(For all data)
# Remember to add the non-overlap data
features = cbind(abs(t(z.decor[2:10,])))

for (i in 1:length(files_no)){features = cbind(features, files_no[[i]]$z)}
z.decor.target.z = cbind(z.decor[1,])

colnames(features) = c('adhd2019','alco2018','dx.dep','ed2019','mdd2019','ocd2018','ptsd2019','tense','wal.scz2018','ad','asd2019','cannabis','college2013','edu','ever.csh','ever.sh','ever.st','life.threat.accident','mddco','neuroti2018','no.dep','physical.crime','psycho.experi','reac','scz.eas2019','sleep.dura','wit.death')

# (6) lasso for chosing the best variables
lasso.cv.bio = cv.glmnet(features, abs(z.decor.target.z), family='gaussian', nlambda=50, alpha=1, standardize=T, parallel=T)
coef_list.min.bio = as.matrix(coef(lasso.cv.bio, s = "lambda.min"))
coef_list.min.bio = unlist(coef_list.min.bio)[-1] 
z.combine.lasso = features[, coef_list.min.bio != 0]
	
cat('Variables not selected from GWAS:', setdiff(colnames(features), colnames(z.combine.lasso)), '\n')

fdr.st = FDRreg(z.decor.target.z, z.combine.lasso, nulltype = 'theoretical', method = 'pr')
fdr.st2 = FDRreg(z.decor.target.z, z.combine.lasso, nulltype = 'empirical', method = 'pr')

# (7) Assessment-Features
library(HelpersMG)
features.se = SEfromHessian(fdr.st$model$hessian)
features.coef = fdr.st$model$coef
features.z.score = features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
features.pvalue = 2*pnorm(abs(features.z.score), lower.tail = F)
assessment = cbind(colnames(features), features.pvalue, features.coef[c(2:length(features.coef))], features.se[c(2:length(features.se))])

# (8) Get the qval
target$qval = p.adjust(target$pval, method="fdr", n=length(target$pval))
target$FDR.the = fdr.st$FDR
target$FDR.emp = fdr.st2$FDR
fwrite(as.data.frame(target), 'bd2018.fdrreg.txt',sep=',')

# (9) Compare FDR-reg result
thresholds = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.001, 5e-04, 5e-06, 5e-08)

count_below_threshold = function(data, thresholds) {
	sapply(thresholds, function(th) length(which(data < th)))
}

qval = count_below_threshold(target$qval, thresholds)
fdr.the = count_below_threshold(target$FDR.the, thresholds)
fdr.emp = count_below_threshold(target$FDR.emp, thresholds)

sum_data = data.frame(qval, fdr.the, fdr.emp)
rownames(sum_data) = c('< 0.5', '< 0.4', '< 0.3', '< 0.2', '< 0.1', '< 0.05', '< 0.04', '< 0.03', '< 0.02', '< 0.01', '< 0.001', '< 5e-04', '< 5e-06', '< 5e-08')

# (9) Prepare data for MAGMA
setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018')

# Target
target.decor.z = cbind(z.decor[1,], target)
colnames(target.decor.z)[1] = 'z.decor'
target.decor.z$p.decor = 2 * pnorm(abs(target.decor.z$z.decor), lower.tail = F)
fwrite(as.data.frame(target.decor.z), 'bd2018.overlap.4magma.txt', sep = ',')

# Variants with sample overlap
for (i in c(1:length(files_with))){
	with.decor.z = cbind(z.decor[i+1,], files_with[[i]])
	colnames(with.decor.z)[1] = 'z.decor'
	name = listnames_with[i]
	with.decor.z$p.decor = 2 * pnorm(abs(with.decor.z$z.decor), lower.tail = F)
    fwrite(as.data.frame(with.decor.z), paste0(name,'.overlap.4magma.txt'), sep=',')
	}

# Variants with no sample overlap
	for (i in c(1:length(files_no))){
	no.decor.z = files_no[[i]]
	name = listnames_no[i]
	no.decor.z$p.decor = no.decor.z$pval
	no.decor.z$z.decor = no.decor.z$z
    fwrite(as.data.frame(no.decor.z), paste0(name, '.overlap.4magma.txt'), sep=',')
	}

#--- quit and save result ----#
fwrite(as.data.frame(assessment), '/exeh_4/jinghong_qiu/SO_Lab/Contribution/snp-bd2018.csv', sep=',')
fwrite(as.data.frame(sum_data), '/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/snp_bd2018.csv', sep=',')
save(fdr.st, fdr.st2, file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/snp_bd2018.RData")

# Time track
end_time = Sys.time()
cat('Time Start:', format(start_time, "%a %b %d %X %Y"), '\n')
cat('Time End:', format(end_time, "%a %b %d %X %Y"), '\n')
cat('Time consuming:', end_time - start_time, '\n')

q('no')
