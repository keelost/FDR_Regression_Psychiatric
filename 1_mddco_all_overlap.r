#Snp overlapping#
	
	#R-353
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_2_lasso/mddco')
	library(data.table);library(stringr);library(dplyr)	
	
#--------------Loading data------------#	
	target = fread('/exeh_4/jinghong_qiu/SO_Lab/impute/mddco.impute.map.txt')
	listnames = c('bd2018','cannabis','dx.dep','ed2019','edu','inte','neuroti2018','ocd2018','scz.eas2019','wal.scz2018','mdd2019')
	files=lapply(paste0('/exeh_4/jinghong_qiu/SO_Lab/impute/',listnames,'.impute.map.txt'),fread,stringsAsFactor=F)

#--------------Get the list of overlapping snps------------#
	all_snps = target$snpid
	for (rank in 1:length(files)){
	all_snps = intersect(all_snps,files[[rank]]$snpid)
	}
	cat('number of snps:',length(all_snps))
	for (rank in 1:length(files)){
	name = listnames[rank]
	trait = files[[rank]]
	idx = match(all_snps,trait$snpid)
	trait.overlap = trait[idx,]
	fwrite(as.data.frame(trait.overlap),paste0(name,'.overlap.txt'),sep=' ')
	}
	idx2 = match(all_snps,target$snpid)
	target.overlap = target[idx2,]
	fwrite(as.data.frame(target.overlap),'mddco.overlap.txt',sep=' ')
	
	q('no')
	
	 