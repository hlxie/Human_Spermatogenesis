## Introduction 


Research on epigenetic regulation of human spermatogenesis remains limited, particularly in understanding the epigenetic dysregulation associated with abnormal spermatogenesis. In this study, we employed high-quality scATAC-seq to investigate the epigenetic regulatory states of male germ cells during both normal spermatogenesis and its disruption in non-obstructive azoospermia (NOA) patients. We identified key transcription factors (TFs) and their target gene networks, demonstrating that their synchronized and wave-like activities are critical for driving the highly ordered progression of spermatogenesis. Additionally, we uncovered significant epigenetic changes at DNA double-strand break (DSB) hotspots, tracing their formation and repair processes at the single-cell level. We showed that a male germ cell set up open chromatin state for several thousands of DSB hotspots. Of these, several hundreds go on to forms DSBs and finally several dozens of these DSBs forms crossovers, permitting proper exchange of genetic materials between parental alleles through homologous recombination. In NOA patients, spermatogenesis is arrested at the zygotene stage with many master TFs showing aberrant dephasing and decoupling. These regulatory failures hinder the spermatogenesis from progressing to the next developmental stage properly. Additionally, we identified point mutations as genetic basis for the dysregulation of TF motif activation such as for NFYB-PITX2-HMGA1 axis. Abnormal DSB hotspot dynamics were also revealed in NOA patients, with approximately half of the opened DSB hotspots in each cell lost accessibility compared to normal spermatogenesis, contributing to a 2-3 fold increase in aneuploid sperm cells. Our findings provide a deeper view of key gene regulatory mechanisms underlying human spermatogenesis and offers new perspectives on the pathogenesis of NOA.

## Scripts
### ./analysis directory
> *01.preprocessing_pipeline.java*

This is a Nextflow script designed for performing essential quality control and alignment of raw scATAC-seq data, ultimately generating fragments files.

> *02.OA_ArchR_pip.R*

We utilized samples from OA patients to construct a chromatin accessibility atlas for normal spermatogenesis. The fragments files generated in the previous step were used as input for ArchR to perform the following analyses: **quality control**, **cell clustering**, **visualization**, **peak calling**, **GeneScore calculation**, **differential GeneScore analysis**, **differential accessibility peak (DAP) analysis**, **enrichment analysis**, **TF motif activity analysis**, and **trajectory analysis**. These analyses correspond to the relevant content in Figures 1 and 2.

> *03.OA_NOA_ArchR_pip.R*

Similar to the script `02.OA_ArchR_pip.R`, this analysis incorporates cells from NOA samples, and the cell types of NOA samples are inferred and identified based on the labels from OA samples.These analyses correspond to the relevant content in **Figures 4 and 5**.

> *04.Joint_scRNA-seq_and_scATAC-seq_analysis.R*

Integrated analysis and visualization of scATAC-seq and scRNA-seq data were performed on OA samples. These analyses correspond to the relevant content in **Figures 2**.

> *05.Pando_analysis.R*

Pando was used to construct the transcription factor (TF) regulatory network, and the network was pruned as needed to generate the `ALL_network_positive_negative.csv` file. These analyses correspond to the relevant content in **Figures 2**.

> *06.DSB_hotspots.R*

Human DMC1 ChIP data was overlapped with chromatin accessible regions of L and Z-stage spermatocytes to generate the `LZ_DMC1.bed` file. The cell-peak binary matrix was extracted, and peaks and cell IDs to be retained were filtered for heatmap visualization. DNA double-strand break (DSB) hotspots were grouped into four clusters, followed by enrichment analysis. These analyses correspond to the relevant content in **Figures 3**.

> *07.Identification_of_sperm_aneuploidy.R*

CNV (copy number variation) was calculated for each cell type using a 5Mb window with OA samples as controls, excluding the hg19 blacklist regions. Each cell type had its own corresponding control group. These analyses correspond to the relevant content in **Figures 6**.

> *08.NOA_call_snv.R*

The OA data was used to construct the control file `OA_pon.vcf.gz`. Mutect2 was then used to calculate SNVs (single nucleotide variants) for single cells and merged SNVs from all cells of each patient. Finally, only SNVs supported by three or more cells and with a population frequency below 0.001 were retained. These analyses correspond to the relevant content in **Figures 7**.

> *09.LD_score_regression.sh*

Statistical significance and linkage disequilibrium information from GWAS data were used to infer the association between genetic variants and cell types.

> *10.Chromafold.sh*

Chromafold was employed to predict the Hi-C structure of normal spermatocytes during the Z stage. These analyses correspond to the relevant content in **Figures 7**.




