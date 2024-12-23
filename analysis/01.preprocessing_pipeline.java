#!/usr/bin/env nextflow

ref= "Path to hg19.fa"
cutadapt= "Path to cutadapt software"
bowtie2= "Path to bowtie2 software"
sinto= "Path to sinto software"


params.raw_dir= "Path to single cell raw data"
params.clean_dir= "01.Clean_data"
params.bam_dir= "02.Alignment"
params.frag_dir= "03.Single_frag"
params.merge_frag_dir= "04.Merge_frag"


params.singlecell_raw="${params.raw_dir}/*fq.gz"
Channel
        .fromPath(params.singlecell_raw)
        .ifEmpty { error "Cannot find any reads matching" }
        .set {each_raw}
        .println()

process trim{

        label "cnlong"
        publishDir params.clean_dir, mode: 'copy', overwrite: false

        input:
                each file(scRaw_cell) from each_raw

        output:
                file("*trimed.*fq.gz") into scClean_data
                file("*.txt") into trim_result

        script:
                """
                ${cutadapt} \
                    -q 20 \
                        -g AGATGTGTATAAGAGACAG \
                        -G AGATGTGTATAAGAGACAG \
                        -a CTGTCTCTTATACACATCT \
                        -A CTGTCTCTTATACACATCT \
                        --times 1 \
                        --minimum-length 30 \
                        -o ${scRaw_cell.getSimpleName()}.trimed.R1.fq.gz \
                        -p ${scRaw_cell.getSimpleName()}.trimed.R2.fq.gz \
                        ${scRaw_cell.getSimpleName()}.R1.fq.gz \
                        ${scRaw_cell.getSimpleName()}.R2.fq.gz \
                        > ./log/log.txt

                """
 }



process align{
    label "cnlong"
    publishDir params.bam_dir, mode: 'copy', overwrite: true

        input:
            file(scClean_cell) from scClean_data

        output:
                file("*.archr.sorted.bam") into Archr_bam
                file("*.archr.sorted.bam.bai") into Archr_bai

        script:
                """
                ${bowtie2} \
                        --mm --no-unal \
                        --very-sensitive \
                        -x ${ref} \
                        -1 ${scClean_cell.getSimpleName()}.trimed.R1.fq.gz \
                        -2 ${scClean_cell.getSimpleName()}.trimed.R2.fq.gz \
                        | samtools view -bS -q 30 \
                        > ${scClean_cell.getSimpleName()}.mapQ30.bam


                samtools view ${scClean_cell.getSimpleName()}.mapQ30.bam -h \
                        | awk -vOFS='\t' '{ print $0,"CB:Z:"${scClean_cell.getSimpleName()} }' sample=${scClean_cell.getSimpleName()} \
                        | samtools view -bS - \
                        > ${scClean_cell.getSimpleName()}.archr.bam

                samtools sort \
                        -o ${scClean_cell.getSimpleName()}.archr.sorted.bam \
                        ${scClean_cell.getSimpleName()}.archr.bam


                samtools index \
                    ${scClean_cell.getSimpleName()}.archr.sorted.bam


                """

}



process single_fragments{
    label "cnlong"
    publishDir params.frag_dir, mode: 'copy', overwrite: true

    input:
            file(scArchr_bam) from Archr_bam
            file(scArchr_bai) from Archr_bai

        output:
                file("*sinto.fragments.bed") into fragments_bed

        script:
                """
                ${sinto} fragments \
                        -b ${scArchr_bam.getSimpleName()}.archr.sorted.bam \
                    -f ${scArchr_bam.getSimpleName()}.sinto.fragments.bed

                 """

}

process merge_batch_fragments{
     label "cnlong"
         publishDir params.merge_frag_dir, mode: 'copy', overwrite: true

         input:
         file(ac_fragments) from fragments_bed.collect()

         output:
             file("batch.fragments.sorted.bed.gz") into All_fragments_bed
                 file("batch.fragments.sorted.bed.gz.tbi") into All_fragments_bed_index

        script:
           """
           cat *sinto.fragments.bed \
               > batch.fragments.bed

           sort -k 1,1 -k2,2n \
                   batch.fragments.bed \
                   > batch.fragments.sorted.bed

           bgzip batch.fragments.sorted.bed
           tabix batch.fragments.sorted.bed.gz



       """

}

