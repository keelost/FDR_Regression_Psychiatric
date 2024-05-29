#Snp overlapping#
	
	#R-353
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/bd2012')
	library(data.table);library(stringr);library(dplyr)	
	
#--------------Loading data------------#	
	target = fread('/exeh_4/jinghong_qiu/SO_Lab/impute/Special_Targets/bd2012.impute.map.txt')
	listnames = c('ad','alco2018','asd2019','cannabis','college2013','dx.dep','ed2019','edu','ever.csh','ever.sh','ever.st','life.threat.accident','mdd2019','mddco','neuroti2018','no.dep','ocd2018','ptsd2019','scz.eas2019','sleep.dura','wal.scz2018')
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
	fwrite(as.data.frame(target.overlap),'bd2012.overlap.txt',sep=' ')
	
	q('no')
	
	 