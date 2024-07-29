library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(ArchR)
library(GenomicRanges)
library(ChIPseeker)
phastConsElements20Mammals.UCSC.hg19 <- readPeakFile("/date/xiehaoling/software/Pando/Pando-main/data/phastConsElements20Mammals.UCSC.hg19.bed")
library(tidyverse)
library(Pando)




pbmc_rna <- readRDS("/date/xiehaoling/Sperm_atac/OA_Y/RNA/OA_RNA.rds")
pbmc <- readRDS("/date/xiehaoling/Sperm_atac/OA_Y/Pando/Pando_peak.rds")

DefaultAssay(pbmc) <- 'Atac_RNA'
transfer.anchors <- FindTransferAnchors(
                                                reference = pbmc_rna,
                                                query = pbmc,
                                                features = unique(c(VariableFeatures(object = pbmc_rna))),
                                                #reference.assay = "RNA", query.assay = "peaks",
                                                reduction = 'cca'
                                                        )


pbmc_rna=FindVariableFeatures(pbmc_rna, selection.method = "vst", nfeatures = 4000)
genes.use <- VariableFeatures(pbmc_rna)
refdata <- GetAssayData(pbmc_rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc[["lsi"]],dims =2:30)
pbmc[["RNA"]] <- imputation
pbmc[['RNA']]=FindVariableFeatures(pbmc[['RNA']], selection.method = "vst", nfeatures = 4000)


saveRDS(pbmc, file = "/date/xiehaoling/Sperm_atac/OA_Y/Pando/Pando_need_RNA_ATAC.rds")

###  pando  ###

DefaultAssay(pbmc) <- 'peaks'

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

muo_data <- initiate_grn(
                                                 pbmc,
                                                 rna_assay = 'RNA',
                                                 peak_assay = 'peaks',
                                                 regions = phastConsElements20Mammals.UCSC.hg19

                                                 )

load('/date/xiehaoling/software/Pando/Pando-main/data/motifs.RData')
load('/date/xiehaoling/software/Pando/Pando-main/data/motif2tf.RData')

library(BSgenome.Hsapiens.UCSC.hg19)
muo_data <- find_motifs(
                                                muo_data,
                                                pfm = motifs,
                                                motif_tfs = motif2tf,
                                                genome = BSgenome.Hsapiens.UCSC.hg19
                                                )

regions <- NetworkRegions(muo_data)

library(doParallel)
#registerDoParallel(4)
muo_data <- infer_grn(
                                        muo_data,
                                        genes = genes.use,
                                        peak_to_gene_method = 'Signac',
                                        parallel = F
                                         )


GetNetwork(muo_data)

muo_data2 <- find_modules(
                                        muo_data,
                                        p_thresh = 0.01,
                                        nvar_thresh = 10,
                                        min_genes_per_module = 3,
                                        rsq_thresh = 0.1
                                                )

modules <- NetworkModules(muo_data2)
muo_data_gg <- get_network_graph(muo_data2,
                                        features =unique(modules@meta$tf),#unique(unique(TF_name)),
                                        #features =unique(unique(TF_name)),
                                        random_seed = 2188
                                        # umap_method = c("corr"),
                                        )

plot_network_graph(muo_data_gg)

data_network=as.data.frame(coef(muo_data2))

data_network_0.5_0.05=data_network[data_network$padj<0.5 & abs(data_network$corr)>0.05,]

write.table(data_network_0.5_0.05,"/date/xiehaoling/Sperm_atac/OA_Y/Pando/network1/ALL_network_positive_negative.csv",quote = F,row.names = F)

