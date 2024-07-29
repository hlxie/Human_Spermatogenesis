cat $1 | while read line
do
        bed=`(echo $line|awk '{print $1}')`
        mkdir ${bed}_folder  #${bed}.bed is the cell type-specific peak region.
        for i in {1..22}
        do
                python2 /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/make_annot.py \
                    --bed-file ${bed}.bed \
                    --bimfile /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/EUR/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim     \
                    --annot-file ${bed}_folder/chr${i}.annot.gz
                gunzip -c ${bed}_folder/chr${i}.annot.gz |awk 'BEGIN{OFS="\t"}{if(NR==1){print "ALL",$1}else{print 1,$1}}' |gzip -c > ${bed}_folder/Add_all.chr${i}.annot.gz
        done
done





cat $1 | while read line
do
        bed=`(echo $line|awk '{print $1}')`
        for i in {1..22}
        do
                    python2 /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/ldsc.py \
                    --l2 \
                    --bfile /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/EUR/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} \
                    --ld-wind-cm 0.1 \
                    --annot ${bed}_folder/Add_all.chr${i}.annot.gz \
                    --thin-annot \
                    --out  ${bed}_folder/Add_all.chr${i}
        done
done

python2 /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/munge_sumstats.py \
        --sumstats Hg19_GCST90239718_build \
        --merge-alleles /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/w_hm3.snplist \
        --chunksize 500000 \
        --out NOA \
            --N 2500 \
            --frq Q

python2 /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/munge_sumstats.py \
        --sumstats Hg19_GCST90239717_build \
        --merge-alleles /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/w_hm3.snplist \
                --chunksize 500000 \
                --out NOA \
                --N 2500 \
                --frq Q

python2 /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/munge_sumstats.py \
        --sumstats Hg19_GCST90239716_build \
        --merge-alleles /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/w_hm3.snplist \
        --chunksize 500000 \
        --out NOA \
        --N 2500 \
        --frq Q

cat $1 | while read line
do
 bed=`(echo $line|awk '{print $1}')`
 python2  /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/ldsc.py \
        --h2 Hg19_GCST90239716_NOA.sumstats.gz \
        --ref-ld-chr ${bed}_folder/Add_all.chr \
        --overlap-annot \
        --n-blocks 1000 \
            --print-coefficients \
            --w-ld-chr /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
            --frqfile-chr /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/1000G_Phase3_frq/1000G.EUR.QC. \
            --out GCST90239716_NOA_${bed}

 python2  /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/ldsc.py \
                --h2 Hg19_GCST90239717_NOA.sumstats.gz \
                --ref-ld-chr ${bed}_folder/Add_all.chr \
                --overlap-annot \
                --n-blocks 1000 \
                --print-coefficients \
                --w-ld-chr /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                --frqfile-chr /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/1000G_Phase3_frq/1000G.EUR.QC. \
                --out GCST90239717_NOA_${bed}

 python2  /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/LDSC/ldsc-master/ldsc.py \
                --h2 Hg19_GCST90239718_NOA.sumstats.gz \
                --ref-ld-chr ${bed}_folder/Add_all.chr \
                --overlap-annot \
                --n-blocks 1000 \
                --print-coefficients \
                --w-ld-chr /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                --frqfile-chr /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/Sperm_atac/GWAS/hg19_LD_score/1000G_Phase3_frq/1000G.EUR.QC. \
                --out GCST90239718_NOA_${bed}

done





