#---------------------------------------------#
#---- FDR-regression/SMultiXcan for BD2018----#
#---------------------------------------------#

# SoftWare: R-4.3.3
# Platform: Server 203
# Updated: 2024/08/02
# Author: Jinghong QIU/Shitao RAO/Hon-Cheong SO

#---- load packages and set dir ----#
rm(list = ls())
chooseCRANmirror(ind=87) # for HK mirror setting, getCRANmirrors() can check the number
options(stringsAsFactors = F)
start_time = Sys.time() # time tracking
set.seed(102) # for repeated result
# load packages
library(FDRreg); library(data.table); library(powerplus); library(stringr); library(dplyr); library(glmnet); library(HelpersMG)

#---- set dir ----#
setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2018')

#---- Pipeline ----#
# (1) Get biological information from online resource and merge with gene list
# DAVID Bioinformatics
scz.gene.david = fread(file ='/exeh_4/jinghong_qiu/SO_Lab/magma-library-uniq-ensembl.csv')

# (2) Load data
# Read all predixcan results in a region
filenames = list.files(path='.', pattern = 'allbrain', full.names = T)
files = lapply(filenames, fread, stringsAsFactor=F)

# (3) Analysis
target_rank = which(str_detect(filenames, 'bd2018'))

# Delete NA in pvalue column, generate zscore from smultixcan pvalue and order datasets with gene ID
for (i in c(1:length(files))) {
	files[[i]] = files[[i]][complete.cases(files[[i]]$pvalue), ]
	files[[i]] = files[[i]][order(files[[i]]$gene), ]
	files[[i]]$zscore = qnorm(files[[i]]$pvalue/2)
}
other_file = c(2:length(files))
other_file = other_file[-(target_rank-1)]

# Delete NA in pvalue column, generate zscore from smultixcan pvalue and order datasets with gene ID
if (!identical(files[[target_rank]]$gene, files[[length(files)]]$gene)){
	idx = match(files[[target_rank]]$gene,files[[length(files)]]$gene)
	omit_rank = setdiff(seq(1:tail(idx, n=1)), idx)
	files[[1]] = files[[1]][-omit_rank, ]

	for (i in other_file){
		files[[i]] = files[[i]][-omit_rank, ]
	}
}
	
# Extract target zscore and construct variables z.score matrix
z.combine = abs(files[[1]]$zscore)
for (i in other_file){
	z.combine = cbind(z.combine, abs(files[[i]]$zscore))
}
colnames(z.combine)=c('ad', 'adhd2019', 'alco2018', 'asd2019', 'cannabis', 'college2013', 'dx.dep', 'ed2019', 'edu', 'ever.csh', 'ever.sh', 'ever.st', 'life.threat.accident', 'mdd2019', 'mddco', 'neuroti2018', 'no.dep', 'ocd2018', 'physical.crime', 'psycho.experi', 'ptsd2019', 'reac', 'scz.eas2019', 'sleep.dura', 'tense', 'wal.scz2018', 'wit.death')
target.z = abs(qnorm(files[[target_rank]]$pvalue/2)) * sign(rnorm(length(files[[target_rank]]$pvalue)))

# Generate geneID and output
gene.list = files[[target_rank]]$gene
gene.list2 = strsplit(gene.list, '[.]')
listdata = vector()
for (line in 1:length(gene.list2)){
	listdata[line] = gene.list2[[line]][1]
}

# (3) Perform FDRreg with two different models
fdr.theore = FDRreg(target.z, z.combine, nulltype = 'theoretical', method = 'pr')
fdr.empiri = FDRreg(target.z, z.combine, nulltype = 'empirical', method = 'pr')

# (4) Accessment
features.se = SEfromHessian(fdr.theore$model$hessian)
features.coef = fdr.theore$model$coef
features.z.score = features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
features.pvalue = 2*pnorm(abs(features.z.score), lower.tail=F)
assessment = data.frame(colnames(z.combine), features.pvalue, features.coef[c(2:length(features.coef))], features.se[c(2:length(features.se))])
colnames(assessment) = c('feature', 'p1', 'beta1', 'se1')

# (5) Biological info merge
target.gene.order = files[[target_rank]]
target.gene.order$gene = listdata
colnames(target.gene.order) = gsub("gene", "ENSEMBL_GENE_ID", colnames(target.gene.order))
scz.bio = left_join(target.gene.order[, 1:2], scz.gene.david, by='ENSEMBL_GENE_ID')
scz.bio[is.na(scz.bio)] = 0
bio.use = as.matrix(scz.bio[, c(3:23)])
z.combine.bio=as.matrix(cbind(abs(z.combine), abs(bio.use)))

fdr.target.theore = FDRreg(target.z, z.combine.bio, nulltype = 'theoretical', method = 'pr')
fdr.target.empiri = FDRreg(target.z, z.combine.bio, nulltype = 'empirical', method = 'pr')

# Assessment-fdr.reg-bio
	features.se = SEfromHessian(fdr.target.theore$model$hessian)
	features.coef = fdr.target.theore$model$coef
	features.z.score = features.coef[c(2:length(features.coef))]/features.se[c(2:length(features.se))]
	features.pvalue = 2 * pnorm(abs(features.z.score), lower.tail=F)
	assessment2 = data.frame(colnames(z.combine.bio), features.pvalue, features.coef[c(2:length(features.coef))], features.se[c(2:length(features.se))])
	colnames(assessment2) = c('feature', 'p2', 'beta2', 'se2')
	assessment_sum = merge(assessment, assessment2, by='feature', all=T)

# Combind results
target.gene.order$qval = p.adjust(files[[target_rank]]$pvalue, method='fdr', n=length(files[[target_rank]]$pvalue))
target.gene.order$FDR.the = fdr.theore$FDR
target.gene.order$FDR.emp = fdr.empiri$FDR
target.gene.order$bio.FDR.the = fdr.target.theore$FDR
target.gene.order$bio.FDR.emp = fdr.target.empiri$FDR

thresholds = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.001, 5e-04, 5e-06, 5e-08)

count_below_threshold = function(data, thresholds) {
	sapply(thresholds, function(th) length(which(data < th)))
}

qval = count_below_threshold(target.gene.order$qval, thresholds)
fdr.the = count_below_threshold(target.gene.order$FDR.the, thresholds)
fdr.emp = count_below_threshold(target.gene.order$FDR.emp, thresholds)
bio.fdr.the = count_below_threshold(target.gene.order$bio.FDR.the, thresholds)
bio.fdr.emp = count_below_threshold(target.gene.order$bio.fdr.emp, thresholds)
	       
sum_data = data.frame(qval, fdr.the, fdr.emp, bio.fdr.the, bio.fdr.emp)
rownames(sum_data) = c('< 0.5', '< 0.4', '< 0.3', '< 0.2', '< 0.1', '< 0.05', '< 0.04', '< 0.03', '< 0.02', '< 0.01', '< 0.001', '< 5e-04', '< 5e-06', '< 5e-08')

cat('----------------------- \n All gene number is: ', length(gene.list), sep = "\n")

#--- quit and save result ----#
fwrite(as.data.frame(target.gene.order), file = 'bd2018.finalresult.txt', sep=',')
fwrite(as.data.frame(assessment_sum), '/exeh_4/jinghong_qiu/SO_Lab/Contribution/smulti-bd2018.csv', sep = ',')
fwrite(as.data.frame(sum_data), '/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/smulti-bd2018.csv', sep = ',')	
save(fdr.theore, fdr.empiri, fdr.target.theore, fdr.target.empiri, file = "/exeh_4/jinghong_qiu/SO_Lab/level2/FDR_models/smulti-bd2018.RData")

# Time track
end_time = Sys.time()
cat('Time Start:', format(start_time, "%a %b %d %X %Y"), '\n')
cat('Time End:', format(end_time, "%a %b %d %X %Y"), '\n')
cat('Time consuming:', end_time - start_time, '\n')

q('no')
