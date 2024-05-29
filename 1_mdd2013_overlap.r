#Snp overlapping#
	
	#R-353
	setwd('/exeh_4/jinghong_qiu/SO_Lab/level2/overlapping_1_all/mdd2013')
	library(data.table);library(stringr);library(dplyr)	
	
#--------------Loading data------------#	
	target = fread('/exeh_4/jinghong_qiu/SO_Lab/impute/Special_Targets/mdd2013.impute.map.txt')
	listnames = c('adhd2019','alco2018','antisocial','anx','asd2019','bd2018','cannabis','chrono','dx.dep','ed2019','ever.csh','ever.sh','ever.st','insom2019','inte','life.threat.accident','neo.c','neuroti2018','no.dep','physical.crime','psycho.experi','ptsd2019','scz.eas2019','st','st.dep','wal.scz2018')
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
	fwrite(as.data.frame(target.overlap),'mdd2013.overlap.txt',sep=' ')
	
	q('no')
	
	 