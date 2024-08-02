#----------------------------------------#
#---- FDR-regression/MAGMA for BD2018----#
#----------------------------------------#

# SoftWare: R-4.3.3
# Platform: Server 203
# Updated: 2024/08/02
# Author: Jinghong QIU/Shitao RAO/Hon-Cheong SO

#---- load packages and set dir ----#
rm(list = ls())
chooseCRANmirror(ind=25) # for HK mirror setting, getCRANmirrors() can check the number
options(stringsAsFactors = F)
start_time = Sys.time() # time tracking
set.seed(25) # for repeated result
# load packages
library(FDRreg);library(data.table);library(powerplus);library(stringr);library(dplyr);library(glmnet);library(HelpersMG)

#---- set dir ----#
setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result')

#---- Pipeline ----#
# (1) Load data
# Target
target.gene = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result/bd2018.genes.out')
target.gene.order = target.gene[order(target.gene$GENE),]
# Variables
listnames = c('ad', 'adhd2019', 'alco2018', 'asd2019', 'cannabis', 'college2013', 'dx.dep', 'ed2019', 'edu', 'ever.csh', 'ever.sh', 'ever.st', 'life.threat.accident', 'mdd2019', 'mddco', 'neuroti2018', 'no.dep', 'ocd2018', 'physical.crime', 'psycho.experi', 'ptsd2019', 'reac', 'scz.eas2019', 'sleep.dura', 'tense', 'wal.scz2018', 'wit.death')
files = lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result/', listnames, '.genes.out'), fread, stringsAsFactor = F)

# (2) Data preproccess
fileout.order = files[[1]][order(files[[1]]$GENE),]
z.combine = abs(qnorm(fileout.order$P/2))
for (i in 2:length(files)){
	geneout = files[[i]]
	geneout.order = geneout[order(geneout$GENE),]
	z.combine = cbind(z.combine,abs(qnorm(geneout.order$P/2)))
	}
colnames(z.combine) = listnames

# Extract target z.stat
target.z = abs(qnorm(target.gene.order$P/2)) * sign(rnorm(length(target.gene.order$P)))

# Generate geneID and output for biological annotation
gene.list = target.gene.order$GENE

# (3) Perform FDR regression of 2 nulltypes
fdr.target.theore = FDRreg(target.z, z.combine, nulltype = 'theoretical', method = 'pr')
fdr.target.empiri = FDRreg(target.z, z.combine, nulltype = 'empirical', method = 'pr')

# (4) Assessment
features.se = SEfromHessian(fdr.target.theore$model$hessian)
features.coef = fdr.target.theore$model$coef
features.z.score = features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
features.pvalue = 2 * pnorm(abs(features.z.score), lower.tail=F)
assessment = data.frame(colnames(z.combine), features.pvalue, features.coef[c(2:length(features.coef))], features.se[c(2:length(features.se))])
colnames(assessment) = c('feature', 'p1', 'beta1', 'se1')

# (5) get biological information from online resource and merge with gene list
# DAVID Bioinformatics
scz.gene.david = fread('/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-entrez.csv')
colnames(target.gene.order) = gsub("GENE", "ID", colnames(target.gene.order))
scz.bio = left_join(target.gene.order[, 1:3], scz.gene.david, by='ID')
scz.bio[is.na(scz.bio)] = 0
scz.bio.use = scz.bio[, c(5:25)]
z.combine.bio = as.matrix(cbind(abs(z.combine), abs(scz.bio[, c(5:25)])))
fdr.target.theore.bio = FDRreg(target.z, z.combine.bio, nulltype = 'theoretical', method = 'pr')
fdr.target.empiri.bio = FDRreg(target.z, z.combine.bio, nulltype = 'empirical', method = 'pr')

# (6) Combind results
target.gene.order$FDR.the = fdr.target.theore$FDR
target.gene.order$FDR.emp = fdr.target.empiri$FDR
target.gene.order$qval = p.adjust(target.gene.order$P, method="fdr", n=length(target.gene.order$P))
target.gene.order$bio.FDR.the = fdr.target.theore.bio$FDR
target.gene.order$bio.FDR.emp = fdr.target.empiri.bio$FDR

thresholds = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.001, 5e-04, 5e-06, 5e-08)

count_below_threshold = function(data, thresholds) {
	sapply(thresholds, function(th) length(which(data < th)))
}

qval = count_below_threshold(target.gene.order$qval, thresholds)
fdr.the = count_below_threshold(target.gene.order$FDR.the, thresholds)
fdr.emp = count_below_threshold(target.gene.order$FDR.emp, thresholds)
bio.fdr.the = count_below_threshold(target.gene.order$bio.FDR.the, thresholds)
bio.fdr.emp = count_below_threshold(target.gene.order$bio.fdr.emp, thresholds)

cat('----------------------- \n All gene number is: ', length(gene.list), sep = "\n")

#--- quit and save result ----#
fwrite(as.data.frame(gene.list), 'bd2018.gene.list.txt')
fwrite(as.data.frame(target.gene.order), 'bd2018.gene.bio.fdrreg.txt', sep=',')
fwrite(as.data.frame(assessment_sum), '/exeh_4/jinghong_qiu/SO_Lab/Contribution/magma-bd2018.csv', sep=',')
sum_data = data.frame(qval, fdr.the, fdr.emp, bio.fdr.the, bio.fdr.emp)
fwrite(as.data.frame(sum_data), '/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/magma_bd2018.csv', sep=',')	
save(fdr.target.theore, fdr.target.empiri, fdr.target.theore.bio, fdr.target.empiri.bio, file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/magma_bd2018.RData")

# Time track
end_time = Sys.time()
cat('Time Start:', format(start_time, "%a %b %d %X %Y"), '\n')
cat('Time End:', format(end_time, "%a %b %d %X %Y"), '\n')
cat('Time consuming:', end_time - start_time, '\n')

q('no')
