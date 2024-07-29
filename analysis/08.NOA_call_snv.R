
######  OA data is used to construct the control file OA_pon.vcf.gz  ######

cat $1 | while read line
do
sample=`(echo $line|awk '{print $1}')`

java -jar -Xmx4g \
    -Djava.io.tmpdir=./tmp \
    /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/picard.jar \
    MarkDuplicates \
    I=${sample}.bam \
    O=${sample}.markdup.bam \
    M=${sample}_dup.metrics \
    CREATE_INDEX=true \
    ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        TMP_DIR=./tmp

        java -jar -Xmx4g /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/picard.jar BuildBamIndex I=${sample}.markdup.bam

        samtools addreplacerg -r "ID:${bam.getSimpleName()}" -r "LB:${bam.getSimpleName()}" -r "SM:${bam.getSimpleName()}" -r "PL:ILLUMINA" -o ${sample}.MarkDup.RG.bam ${sample}.markdup.bam

        samtools index ${sample}.MarkDup.RG.bam

gatk Mutect2 \
        -R /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/refrence/gatk_hg19/HG19/hg19.fa \
        -I ${sample}.MarkDup.RG.bam \
        -tumor ${sample} \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        -O ${sample}.vcf.gz


done

gatk --java-options "-Xmx4g -Xms4g"  GenomicsDBImport \
                -R /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/refrence/bwa_hg19/hg19_24.fa \
                -L /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/pacbioSeq_old/reference/gatk4/bundle/hg19/interval_list/Genome_interval_list-main/wgs_calling_regions.hg19.interval_list \
                --genomicsdb-workspace-path pon_db \
                -V OA1.vcf.gz \
                -V OA2.vcf.gz \
                -V OA3.vcf.gz \
                -V OA4.vcf.gz \
                -V OA5.vcf.gz \
                -V OA6.vcf.gz \
                -V OA7.vcf.gz \
                -V OA8.vcf.gz \
                -V OA9.vcf.gz \
                -V OA10.vcf.gz \
                -V OA11.vcf.gz \
                -V OA12.vcf.gz \
                -V OA13.vcf.gz \
                -V OA14.vcf.gz \
                -V OA15.vcf.gz \
                -V OA16.vcf.gz
gatk CreateSomaticPanelOfNormals \
     -R /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/refrence/bwa_hg19/hg19_24.fa \
     -V gendb://pon_db \
     -O OA_pon.vcf.gz


######   For each NOA single cell or merged patient file, use Mutect2 to call SNVs   ######
######   Using the merged file of all cells from NOA10 as an example   ######

gatk --java-options "-Xmx2g" Mutect2 \
      -R /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/refrence/bwa_hg19/hg19_24.fa \
      -I NOA10.MarkDup.RG.bam \
      -pon OA_pon.vcf.gz \
      --germline-resource /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/pacbioSeq_old/reference/gatk4/bundle/hg19/af-only-gnomad.hg19_filter.vcf.gz \
      --af-of-alleles-not-in-resource 0.001 \
      --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
      -L /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/pacbioSeq_old/reference/gatk4/bundle/hg19/interval_list/Genome_interval_list-main/wgs_calling_regions.hg19.interval_list \
      -O NOA10_somatic_m2.vcf.gz

/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/anaconda3/envs/ldsc/bin/gatk  FilterMutectCalls \
      -R /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/refrence/bwa_hg19/hg19_24.fa \
      -V NOA10_somatic_m2.vcf.gz \
      -O NOA10_somatic.vcf.gz


bcftools filter -i ' FILTER="PASS" ' NOA10_somatic.vcf.gz > NOA10_PASS_somatic.vcf
bgzip NOA10_PASS_somatic.vcf
bcftools index NOA10_PASS_somatic.vcf.gz


bcftools filter -i ' AF>=0.1 ' NOA10_PASS_somatic.vcf.gz > NOA10_PASS_somatic.AF0.01.vcf
less NOA10_PASS_somatic.AF0.01.vcf| awk '{if($0 ~ /^#/) print $0; else if(length($4)==1) print $0 }' | awk '{if($0 ~ /^#/) print $0; else if(length($5)==1) print $0 }' > NOA10_PASS_somatic.AF0.01.SNV.vcf

bgzip NOA10_PASS_somatic.AF0.01.SNV.vcf
bcftools index NOA10_PASS_somatic.AF0.01.SNV.vcf.gz

###### To calculate SNVs supported by more than three cells, it is necessary to merge all single-cell VCF files.  ######
ls *.vcf.gz > sample.list

merge_sample=""

for sampleid in $(cat sample.list)
do

        merge_sample=$merge_sample" "$sampleid

done

bcftools merge $merge_sample -o NOA10_merged.vcf
bgzip NOA10_merged.vcf


bcftools view -i 'count(GT="mis")<767' NOA10_merged.vcf.gz > More_4_NOA10_merged.vcf ##767+3 is total cell in NOA10
bgzip More_4_NOA10_merged.vcf
bcftools index More_4_NOA10_merged.vcf.gz
bcftools isec -p Single_sudoBulk_overlap -c all -Oz More_4_NOA10_merged.vcf.gz NOA10_PASS_somatic.AF0.01.SNV.vcf.gz


######   To filter out potential SNP sites, we include files from three databases: dbSNP_up_0.001.vcf.gz, dbSNP_COMMON.vcf.gz, and gnomad.exomes.r2.0.2.sites.0.001.vcf.gz. ######
bcftools isec -p DP9_dir_dbSNP_up -c all -Oz 0003.vcf.gz /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/pacbioSeq_old/reference/gatk4/bundle/hg19/dbSNP_all/dbSNP_up_0.001.vcf.gz

bcftools isec -p DP9_dir_dbSNP_up/dbSNP_COMMON -c all -Oz DP9_dir_dbSNP_up/0000.vcf.gz /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/pacbioSeq_old/reference/gatk4/bundle/hg19/dbSNP_all/dbSNP_COMMON.vcf.gz


bcftools isec -p DP9_dir_dbSNP_up/dbSNP_COMMON/gnomad -c all -Oz DP9_dir_dbSNP_up/dbSNP_COMMON/0000.vcf.gz /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/pacbioSeq_old/reference/gatk4/bundle/hg19/gnomad_all/gnomad.exomes.r2.0.2.sites.0.001.vcf.gz


######  Annotate the SNVs   ######
perl /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/ANNOVAR/annovar/convert2annovar.pl -format vcf4 DP9_dir_dbSNP_up/dbSNP_COMMON/gnomad/0000.vcf.gz -out ./NOA.SNV.avinput

perl /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/ANNOVAR/annovar/annotate_variation.pl -geneanno -buildver hg19 ./NOA.SNV.avinput /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/xiehaoling/software/ANNOVAR/annovar/humandb/
