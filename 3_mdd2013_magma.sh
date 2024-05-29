cd /exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013
###注意！！有些Target文件不在library里面
for i in *.4magma.txt
do
    outputname=${i%.o*}
	y=$(echo /exeh_3/rstao/dr.so/002.meta/002.clear.data/library/$outputname'.clear.txt')
    position=$(awk 'BEGIN{ FS=" ";}
    NR == 1 { 
       for (i=1; i <= NF; ++i) {
          if ($i == "n") {
          position=i;
         }
        }
      }
    END{ print position ;}' $y)
    lines=$(awk 'NR!=1{print $var}' var="${position}" $y | head -1 )


nohup /exeh_3/rstao/programs/magma/magma --bfile /exeh_4/jinghong_qiu/Program/magma/g1000_eur --pval $i use=snpid,p.decor N=$lines --gene-annot /exeh_3/rstao/programs/magma/g1000.eur.10k.snp2gene.genes.annot  --out /exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013/Result/$outputname >/exeh_4/jinghong_qiu/logfile/mdd2013/$outputname &


done

/exeh_4/jinghong_qiu/Program/magma/magma --bfile /exeh_4/jinghong_qiu/Program/magma/g1000_eur --pval /exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013/mdd2013.overlap.4magma.txt use=snpid,p.decor N=18759 --gene-annot /exeh_3/rstao/programs/magma/g1000.eur.10k.snp2gene.genes.annot  --out /exeh_4/jinghong_qiu/SO_Lab/level2/magma_1_all/mdd2013/Result/mdd2013 &