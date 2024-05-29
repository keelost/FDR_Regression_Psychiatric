#############
#----FDR----#
#############

#R-353
	library(data.table);library(powerplus)
#载入数据:SCZ
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2012/scz2012.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/wal.scz2018/wal.scz2018.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的2个值
	idx2 = which(old_data$qval<0.05)

	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	#进行overlap
	aaa = intersect(new_change$snpid,old_change_fdrr$snpid)
	bbb = intersect(new_change$snpid,old_change_origin$snpid)
	fold_change = length(aaa)/length(bbb)
	all_overlap = intersect(new_data$snpid,old_data$snpid)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	cat(length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,sep='\n')

#载入数据：SCZ2014
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2012/scz2012.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2014/scz2014.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的2个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	#进行overlap
	aaa = intersect(new_change$snpid,old_change_fdrr$snpid)
	bbb = intersect(new_change$snpid,old_change_origin$snpid)
	fold_change = length(aaa)/length(bbb)
	all_overlap = intersect(new_data$snpid,old_data$snpid)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	cat(length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,sep='\n')

#载入数据:BD
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/bd2018/bd2018.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/bd2012/bd2012.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的2个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	#进行overlap
	aaa = intersect(new_change$snpid,old_change_fdrr$snpid)
	bbb = intersect(new_change$snpid,old_change_origin$snpid)
	fold_change = length(aaa)/length(bbb)
	all_overlap = intersect(new_data$snpid,old_data$snpid)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	cat(length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,sep='\n')

	
	
#载入数据:MDD
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/mdd2013/mdd2013.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/mdd2019/mdd2019.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的2个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	#进行overlap
	aaa = intersect(new_change$snpid,old_change_fdrr$snpid)
	bbb = intersect(new_change$snpid,old_change_origin$snpid)
	fold_change = length(aaa)/length(bbb)
	all_overlap = intersect(new_data$snpid,old_data$snpid)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	cat(length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	
#R-353
	library(data.table);library(powerplus)
#载入数据:SCZ:BIG DATA
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/scz2014/scz2014.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/fdr_1_all/wal.scz2018/wal.scz2018.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的2个值
	idx2 = which(old_data$qval<0.05)
	
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	#进行overlap
	aaa = intersect(new_change$snpid,old_change_fdrr$snpid)
	bbb = intersect(new_change$snpid,old_change_origin$snpid)
	fold_change = length(aaa)/length(bbb)
	all_overlap = intersect(new_data$snpid,old_data$snpid)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	cat(length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,sep='\n')

###############
#----MAGMA----#
###############

#R-353
	library(data.table);library(powerplus)
#载入数据:SCZ
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2012/Result/scz2012.gene.bio.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/wal.scz2018/Result/wal.scz2018.gene.bio.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ID,old_change_fdrr$ID)
	bbb = intersect(new_change$ID,old_change_origin$ID)
	ccc = intersect(new_change$ID,old_change_bio$ID)
	ddd = intersect(new_change$ID,old_change_lasso$ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ID,old_data$ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')
	
#载入数据:SCZ
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2012/Result/scz2012.gene.bio.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2014/Result/scz2014.gene.bio.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ID,old_change_fdrr$ID)
	bbb = intersect(new_change$ID,old_change_origin$ID)
	ccc = intersect(new_change$ID,old_change_bio$ID)
	ddd = intersect(new_change$ID,old_change_lasso$ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ID,old_data$ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')

#载入数据:BD
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2012/Result/bd2012.gene.bio.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/bd2018/Result/bd2018.gene.bio.fdrreg.txt')
	head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ID,old_change_fdrr$ID)
	bbb = intersect(new_change$ID,old_change_origin$ID)
	ccc = intersect(new_change$ID,old_change_bio$ID)
	ddd = intersect(new_change$ID,old_change_lasso$ID)
	all_overlap = intersect(new_data$ID,old_data$ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')
	
#载入数据:MDD
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013/Result/mdd2013.gene.bio.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2019/Result/mdd2019.gene.bio.fdrreg.txt')
	head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ID,old_change_fdrr$ID)
	bbb = intersect(new_change$ID,old_change_origin$ID)
	ccc = intersect(new_change$ID,old_change_bio$ID)
	ddd = intersect(new_change$ID,old_change_lasso$ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ID,old_data$ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(length(all_overlap),length(all_overlap)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(length(all_overlap),length(all_overlap)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(length(all_overlap),length(all_overlap)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')

#载入数据:SCZ:BIG DAT
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/scz2014/Result/scz2014.gene.bio.fdrreg.txt')
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/wal.scz2018/Result/wal.scz2018.gene.bio.fdrreg.txt')
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ID,old_change_fdrr$ID)
	bbb = intersect(new_change$ID,old_change_origin$ID)
	ccc = intersect(new_change$ID,old_change_bio$ID)
	ddd = intersect(new_change$ID,old_change_lasso$ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ID,old_data$ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')
	
	
##################
#----metaxcan----#
##################

#R-353
	library(data.table);library(powerplus)
	brain=c('Brain_Amygdala',
	'Brain_Anterior_cingulate_cortex_BA24',
	'Brain_Caudate_basal_ganglia',
	'Brain_Cerebellar_Hemisphere',
	'Brain_Cerebellum',
	'Brain_Cortex',
	'Brain_Frontal_Cortex_BA9',
	'Brain_Hippocampus',
	'Brain_Hypothalamus',
	'Brain_Nucleus_accumbens_basal_ganglia',
	'Brain_Putamen_basal_ganglia',
	'Brain_Spinal_cord_cervical_c-1',
	'Brain_Substantia_nigra'
	)
	line_all = c('type','fdr.sig','fdrr.sig','fold_change','overlap','big_sig','pvalue','region')
for (brain_rigion in brain){
#载入数据:SCZ(p-value < 2.2e-16  OR 27.07749

	old_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/Result/scz2012.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	new_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/Result/wal.scz2018.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	cat('【',brain_rigion,'】')
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	frame_1 = data.frame(c(length(bbb),length(aaa)),nrow(new_change)-length(bbb),nrow(new_change)-length(aaa))
	frame_2 = data.frame(c(length(bbb),length(ccc)),nrow(new_change)-length(bbb),nrow(new_change)-length(ccc))
	frame_3 = data.frame(c(length(bbb),length(ddd)),nrow(new_change)-length(bbb),nrow(new_change)-length(ddd))
	
	if(all(frame_1>=5))
	{prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop = fisher.test(as.matrix(frame_1),alternative="less")}
	if(all(frame_2>=5))
	{prop2 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop2 = fisher.test(as.matrix(frame_2),alternative="less")}
	if(all(frame_3>=5))
	{prop3 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop3 = fisher.test(as.matrix(frame_3),alternative="less")}
		
	line1 = c('no.bio',length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,brain_rigion)
	line2 = c('bio',length(bbb),length(ccc),fold_change_bio,length(all_overlap),nrow(new_change),prop2$p.value,brain_rigion)
	line3 = c('lasso',length(bbb),length(ddd),fold_change_lasso,length(all_overlap),nrow(new_change),prop3$p.value,brain_rigion)
	line_sum = t(data.frame(line1,line2,line3))
	line_all = rbind(line_all,line_sum)
}
fwrite(as.data.frame(line_all),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/validation-meta-wal.scz2018.csv',sep = ',')
line_all = c('type','fdr.sig','fdrr.sig','fold_change','overlap','big_sig','pvalue','region')
for (brain_rigion in brain){
#载入数据:SCZ2014(p-value < 2.2e-16  OR 27.07749

	old_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/Result/scz2012.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	new_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2014/Result/scz2014.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	cat('【',brain_rigion,'】')
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	
	frame_1 = data.frame(c(length(bbb),length(aaa)),nrow(new_change)-length(bbb),nrow(new_change)-length(aaa))
	frame_2 = data.frame(c(length(bbb),length(ccc)),nrow(new_change)-length(bbb),nrow(new_change)-length(ccc))
	frame_3 = data.frame(c(length(bbb),length(ddd)),nrow(new_change)-length(bbb),nrow(new_change)-length(ddd))
	
	if(all(frame_1>=5))
	{prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop = fisher.test(as.matrix(frame_1),alternative="less")}
	if(all(frame_2>=5))
	{prop2 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop2 = fisher.test(as.matrix(frame_2),alternative="less")}
	if(all(frame_3>=5))
	{prop3 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop3 = fisher.test(as.matrix(frame_3),alternative="less")}
		
	line1 = c('no.bio',length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,brain_rigion)
	line2 = c('bio',length(bbb),length(ccc),fold_change_bio,length(all_overlap),nrow(new_change),prop2$p.value,brain_rigion)
	line3 = c('lasso',length(bbb),length(ddd),fold_change_lasso,length(all_overlap),nrow(new_change),prop3$p.value,brain_rigion)
	line_sum = t(data.frame(line1,line2,line3))
	line_all = rbind(line_all,line_sum)
}
fwrite(as.data.frame(line_all),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/validation-meta-scz2014.csv',sep = ',')
line_all = c('type','fdr.sig','fdrr.sig','fold_change','overlap','big_sig','pvalue','region')

# Brain_Amygdala pvalue 4.962295e-21 OR 13.71524 
# Brain_Anterior_cingulate_cortex_BA24 pvalue 1.754214e-10 OR 28.06867 
# Brain_Caudate_basal_ganglia pvalue 7.846048e-10 OR 34.58714 
# Brain_Cerebellar_Hemisphere pvalue 3.972317e-81 OR 86.36767 
# Brain_Cerebellum pvalue 2.208318e-80 OR 123.5332 
# Brain_Cortex pvalue 1.067127e-84 OR 101.3026 
# Brain_Frontal_Cortex_BA9 pvalue 2.093694e-27 OR 43.87379 
# Brain_Hippocampus pvalue 1.215539e-24 OR 38.59401 
# Brain_Hypothalamus pvalue 3.67169e-32 OR 45.9468 
# Brain_Nucleus_accumbens_basal_ganglia pvalue 8.981974e-56 OR 32.13666 
# Brain_Putamen_basal_ganglia pvalue 1.199119e-54 OR 46.21765 
# Brain_Spinal_cord_cervical_c-1 pvalue 1.767739e-52 OR 34.84072 
# Brain_Substantia_nigra pvalue 9.455935e-37 OR 38.41461 
	
for (brain_rigion in brain){
#载入数据:BD(p-value < 2.2e-16  OR 27.07749

	old_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2012/Result/bd2012.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	new_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2018/Result/bd2018.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	cat('【',brain_rigion,'】')
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
		
	frame_1 = matrix(c(length(bbb),length(aaa),nrow(new_change)-length(bbb),nrow(new_change)-length(aaa)),nrow=2)
	frame_2 = matrix(c(length(bbb),length(ccc),nrow(new_change)-length(bbb),nrow(new_change)-length(ccc)),nrow=2)
	frame_3 = matrix(c(length(bbb),length(ddd),nrow(new_change)-length(bbb),nrow(new_change)-length(ddd)),nrow=2)
	
	if(all(frame_1>=5))
	{prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop = fisher.test(as.matrix(frame_1),alternative="less")}
	if(all(frame_2>=5))
	{prop2 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop2 = fisher.test(as.matrix(frame_2),alternative="less")}
	if(all(frame_3>=5))
	{prop3 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop3 = fisher.test(as.matrix(frame_3),alternative="less")}
		
	line1 = c('no.bio',length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,brain_rigion)
	line2 = c('bio',length(bbb),length(ccc),fold_change_bio,length(all_overlap),nrow(new_change),prop2$p.value,brain_rigion)
	line3 = c('lasso',length(bbb),length(ddd),fold_change_lasso,length(all_overlap),nrow(new_change),prop3$p.value,brain_rigion)
	line_sum = t(data.frame(line1,line2,line3))
	line_all = rbind(line_all,line_sum)
}
fwrite(as.data.frame(line_all),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/validation-meta-bd2018.csv',sep = ',')
line_all = c('type','fdr.sig','fdrr.sig','fold_change','overlap','big_sig','pvalue','region')
# Brain_Amygdala 0.001477087 20.40417 
# Brain_Anterior_cingulate_cortex_BA24 1.811675e-13 8.700587 
# Brain_Caudate_basal_ganglia 1.229594e-06 104.2861 
# Brain_Cerebellar_Hemisphere 1 0 
# Brain_Cerebellum 4.908334e-30 100.4714 
# Brain_Cortex 0.1251772 9.287997 
# Brain_Frontal_Cortex_BA9 0.06283665 19.48306 
# Brain_Hippocampus 1 0 
# Brain_Hypothalamus 1 0 
# Brain_Nucleus_accumbens_basal_ganglia 1.619482e-08 Inf 
# Brain_Putamen_basal_ganglia 1 0 
# Brain_Spinal_cord_cervical_c-1 6.240757e-52 58.07857 
# Brain_Substantia_nigra 3.656417e-10 145.8476 

for (brain_rigion in brain){
#载入数据:MDD(p-value < 2.2e-16  OR 27.07749

	old_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2013/Result/mdd2013.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	new_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2019/Result/mdd2019.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	cat('【',brain_rigion,'】')
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
		
	frame_1 = matrix(c(length(bbb),length(aaa),nrow(new_change)-length(bbb),nrow(new_change)-length(aaa)),nrow=2)
	frame_2 = matrix(c(length(bbb),length(ccc),nrow(new_change)-length(bbb),nrow(new_change)-length(ccc)),nrow=2)
	frame_3 = matrix(c(length(bbb),length(ddd),nrow(new_change)-length(bbb),nrow(new_change)-length(ddd)),nrow=2)
	
	if(all(frame_1>=5))
	{prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop = fisher.test(as.matrix(frame_1),alternative="less")}
	if(all(frame_2>=5))
	{prop2 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop2 = fisher.test(as.matrix(frame_2),alternative="less")}
	if(all(frame_3>=5))
	{prop3 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop3 = fisher.test(as.matrix(frame_3),alternative="less")}
		
	line1 = c('no.bio',length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,brain_rigion)
	line2 = c('bio',length(bbb),length(ccc),fold_change_bio,length(all_overlap),nrow(new_change),prop2$p.value,brain_rigion)
	line3 = c('lasso',length(bbb),length(ddd),fold_change_lasso,length(all_overlap),nrow(new_change),prop3$p.value,brain_rigion)
	line_sum = t(data.frame(line1,line2,line3))
	line_all = rbind(line_all,line_sum)
}
fwrite(as.data.frame(line_all),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/validation-meta-mdd2019.csv',sep = ',')

# Brain_Amygdala 1 0 
# Brain_Anterior_cingulate_cortex_BA24 1 0 
# Brain_Caudate_basal_ganglia 1 0 
# Brain_Cerebellar_Hemisphere 1 0 
# Brain_Cerebellum 1 0 
# Brain_Cortex 1 0 
# Brain_Frontal_Cortex_BA9 0.08641088 3.131522 
# Brain_Hippocampus 1 0 
# Brain_Hypothalamus 1 0 
# Brain_Nucleus_accumbens_basal_ganglia 1 0 
# Brain_Putamen_basal_ganglia 1 0 
# Brain_Spinal_cord_cervical_c-1 1 0 
# Brain_Substantia_nigra 0.1167811 2.708913
line_all = c('type','fdr.sig','fdrr.sig','fold_change','overlap','big_sig','pvalue','region')
	
for (brain_rigion in brain){
#载入数据:SCZ:BIG DATA

	old_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2014/Result/scz2014.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	new_data = read.table(paste0('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/Result/wal.scz2018.finalresult.',brain_rigion,'.test_SetSeed.txt'),sep=' ',header=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	cat('【',brain_rigion,'】')
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
		
	frame_1 = matrix(c(length(bbb),length(aaa),nrow(new_change)-length(bbb),nrow(new_change)-length(aaa)),nrow=2)
	frame_2 = matrix(c(length(bbb),length(ccc),nrow(new_change)-length(bbb),nrow(new_change)-length(ccc)),nrow=2)
	frame_3 = matrix(c(length(bbb),length(ddd),nrow(new_change)-length(bbb),nrow(new_change)-length(ddd)),nrow=2)
	
	if(all(frame_1>=5))
	{prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop = fisher.test(as.matrix(frame_1),alternative="less")}
	if(all(frame_2>=5))
	{prop2 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop2 = fisher.test(as.matrix(frame_2),alternative="less")}
	if(all(frame_3>=5))
	{prop3 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop3 = fisher.test(as.matrix(frame_3),alternative="less")}
		
	line1 = c('no.bio',length(bbb),length(aaa),fold_change,length(all_overlap),nrow(new_change),prop$p.value,brain_rigion)
	line2 = c('bio',length(bbb),length(ccc),fold_change_bio,length(all_overlap),nrow(new_change),prop2$p.value,brain_rigion)
	line3 = c('lasso',length(bbb),length(ddd),fold_change_lasso,length(all_overlap),nrow(new_change),prop3$p.value,brain_rigion)
	line_sum = t(data.frame(line1,line2,line3))
	line_all = rbind(line_all,line_sum)
}
fwrite(as.data.frame(line_all),'/exeh_4/jinghong_qiu/SO_Lab/level2/Result_table/validation-meta-scz.big.csv',sep = ',')




####################
#----SmultiXcan----#
####################

#R-353
	library(data.table);library(powerplus)
#载入数据:SCZ 
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/scz2012.finalresult.txt',fill=T)
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/wal.scz2018.finalresult.txt',fill=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')

#载入数据:SCZ2014
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2012/scz2012.finalresult.txt',fill=T)
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2014/scz2014.finalresult.txt',fill=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')
	
	
#载入数据:BD 47.63656
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2012/bd2012.finalresult.txt',fill=T)
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/bd2018/bd2018.finalresult.txt',fill=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	frame_1 = matrix(c(length(bbb),length(aaa),nrow(new_change)-length(bbb),nrow(new_change)-length(aaa)),nrow=2)
	frame_2 = matrix(c(length(bbb),length(ccc),nrow(new_change)-length(bbb),nrow(new_change)-length(ccc)),nrow=2)
	frame_3 = matrix(c(length(bbb),length(ddd),nrow(new_change)-length(bbb),nrow(new_change)-length(ddd)),nrow=2)
	
	if(all(frame_1>=5))
	{prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop = fisher.test(as.matrix(frame_1),alternative="less")}
	if(all(frame_2>=5))
	{prop2 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop2 = fisher.test(as.matrix(frame_2),alternative="less")}
	if(all(frame_3>=5))
	{prop3 = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')} else
	{prop3 = fisher.test(as.matrix(frame_3),alternative="less")}

	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')

#载入数据:MDD 47.63656
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2013/mdd2013.finalresult.txt',fill=T)
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/mdd2019/mdd2019.finalresult.txt',fill=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')

#载入数据:SCZ 
	old_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/scz2014/scz2014.finalresult.txt',fill=T)
	new_data = fread('/exeh_4/jinghong_qiu/SO_Lab/level2/smultixcan_1_all/wal.scz2018/wal.scz2018.finalresult.txt',fill=T)
	#head(old_data);head(new_data)
	#大数据集sig.fdr
	idx = which(new_data$qval<0.05)
	new_change = new_data[idx,]
	#小数据集的3个值
	idx2 = which(old_data$qval<0.05)
	old_change_origin = old_data[idx2,]
	idx3 = which(old_data$FDR.the<0.05)
	old_change_fdrr = old_data[idx3,]
	idx4 = which(old_data$bio.FDR.the<0.05)
	old_change_bio = old_data[idx4,]
	idx5 = which(old_data$bio.FDR.the.lasso<0.05)
	old_change_lasso = old_data[idx5,]
	#进行overlap
	aaa = intersect(new_change$ENSEMBL_GENE_ID,old_change_fdrr$ENSEMBL_GENE_ID)
	bbb = intersect(new_change$ENSEMBL_GENE_ID,old_change_origin$ENSEMBL_GENE_ID)
	ccc = intersect(new_change$ENSEMBL_GENE_ID,old_change_bio$ENSEMBL_GENE_ID)
	ddd = intersect(new_change$ENSEMBL_GENE_ID,old_change_lasso$ENSEMBL_GENE_ID)
	fold_change = length(aaa)/length(bbb)
	fold_change_bio = length(ccc)/length(bbb)
	fold_change_lasso = length(ddd)/length(bbb)
	all_overlap = intersect(new_data$ENSEMBL_GENE_ID,old_data$ENSEMBL_GENE_ID)
	prop = prop.test(c(length(bbb),length(aaa)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop2 = prop.test(c(length(bbb),length(ccc)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
	prop3 = prop.test(c(length(bbb),length(ddd)),c(nrow(new_change),nrow(new_change)),alternative = 'less')
		
	cat('【no bio】',fold_change,length(bbb),length(aaa),length(all_overlap),nrow(new_change),prop$p.value,sep='\n')
	cat('【bio】',fold_change_bio,length(bbb),length(ccc),length(all_overlap),nrow(new_change),prop2$p.value,sep='\n')
	cat('【lasso】',fold_change_lasso,length(bbb),length(ddd),length(all_overlap),nrow(new_change),prop3$p.value,sep='\n')
