library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(ArchR)


projHeme6=readRDS("/date/xiehaoling/Sperm_atac/OA/Save-OA-ArchR-Project.rds")

use_cell=row.names(projHeme6@cellColData[projHeme6@cellColData$Cell_type %in%c("SSC_1","SSC_2","Diff_ing_SPG","Diff_ed_SPG","PreL","L","Z","PreP","P","D","Second","S1","S2","S3"),])

counts_mat= getMatrixFromProject(projHeme6, useMatrix="PeakMatrix")
PeakMatrix_data= counts_mat@assays@data$PeakMatrix
colnames(PeakMatrix_data)= sapply(strsplit(row.names(projHeme6@cellColData), split = "#"), function(x) x[2])



name_tmp= as.data.frame(ranges(counts_mat@rowRanges))
name_tmp$chr= as.character(seqnames(counts_mat@rowRanges))
row_names= paste(name_tmp$chr, name_tmp$start, name_tmp$end, sep=c(":", "-", ""))
row.names(PeakMatrix_data)= row_names


metadata=as.data.frame(projHeme6@cellColData)

pbmc <- CreateSeuratObject(
                                        counts = PeakMatrix_data,
                                        assay = 'peaks',
                                        project = 'ATAC',
                                        min.cells = 1,
                                        meta.data = metadata
                                                )

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
Annotation(pbmc) <- annotations

pbmc <- RunSVD(pbmc)
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()


gene.activities <- GeneActivity(pbmc)
pbmc[['Atac_RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
                                object = pbmc,
                                assay = 'Atac_RNA',
                                normalization.method = 'LogNormalize',
                                scale.factor = median(pbmc$nCount_Atac_RNA)
                                        )

saveRDS(pbmc, file = "/date/xiehaoling/Sperm_atac/OA_Y/Pando/Pando_peak.rds")

#### Load scRNA-seq data  ###

pbmc_rna <- readRDS("/date/xiehaoling/Sperm_atac/OA_Y/RNA/OA_RNA.rds")
pbmc_rna=FindVariableFeatures(pbmc_rna, selection.method = "vst", nfeatures = 3000)

genes.use <- VariableFeatures(pbmc_rna)
refdata <- GetAssayData(pbmc_rna, assay = "RNA", slot = "data")[genes.use, ]

transfer.anchors <- FindTransferAnchors(reference = pbmc_rna, query = pbmc, features = VariableFeatures(object = pbmc_rna),
                                        reference.assay = "RNA", query.assay = "Atac_RNA", reduction = "cca")


imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc[["lsi"]],
                                                       dims = 2:30)
pbmc[["RNA"]] <- imputation
coembed <- merge(x = pbmc_rna, y = pbmc)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- FindNeighbors(coembed, reduction = "pca", dims = 1:30)
coembed <- FindClusters(coembed, resolution = 1.5)
DimPlot(coembed, group.by = c("orig.ident","seurat_clusters"))


saveRDS(coembed, file = "/date/xiehaoling/Sperm_atac/OA_Y/Pando/coembed.rds")

                                  
