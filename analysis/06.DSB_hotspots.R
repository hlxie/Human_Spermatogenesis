library(ArchR)
library(Seurat)
set.seed(1)

projHeme6=readRDS("/date/xiehaoling/Sperm_atac/OA/Save-OA-ArchR-Project.rds")
LZ_PRDM9_Overlap_peak=read.table("LZ_PRDM9.bed") #The LZ_PRDM9.bed file is an intersection file of DMC1 chip data and the open peaks during the L and Z stages of normal spermatogenesis


PeakMatrix=getMatrixFromProject(projHeme6,useMatrix = "PeakMatrix",binarize = T)
Use_cell=row.names(projHeme6@cellColData[projHeme6@cellColData$Cell_type=="L"| projHeme6@cellColData$Cell_type=="Z",])

PeakMatrix@rowRanges$ID=seq(1:647124)
row_ID=PeakMatrix@rowRanges[PeakMatrix@rowRanges %over% LZ_PRDM9_Overlap_peak]$ID
Use_peak=PeakMatrix@rowRanges[PeakMatrix@rowRanges %over% LZ_PRDM9_Overlap_peak]
Use_peak_name=paste(Use_peak@seqnames,as.data.frame(ranges(Use_peak))$start,as.data.frame(ranges(Use_peak))$end,sep = "_")
PeakMatrix_PRDM9_AA_union=as.matrix(PeakMatrix@assays@data$PeakMatrix[row_ID,Use_cell])
row.names(PeakMatrix_PRDM9_AA_union)=Use_peak_name


DF=as.data.frame(as.matrix(PeakMatrix@assays@data$PeakMatrix[row_ID,row.names(projHeme6@cellColData)]))
row.names(DF)=Use_peak_name
write.table(row.names(DF),"PRDM9_peak.bed",quote = F,row.names = F,sep="\t")



###Remove peaks within 2 kb of the TSS
library(pheatmap)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

PRDM9_peak=readPeakFile("PRDM9_peak.bed")
peakAnno <- annotatePeak(PRDM9_peak, tssRegion=c(-3000, 1000), TxDb=txdb,annoDb="org.Hs.eg.db")
write.csv(as.data.frame(peakAnno@anno[abs(peakAnno@anno$distanceToTSS)>2000,])[,c("seqnames","start","end")],
                            "PRDM9_used_peak_tss_2K.bed",row.names=F,quote=F)



###Trajectory  analysis
Finally_used_peak=read.table("/date/xiehaoling/Sperm_atac/OA_Y/Save-OA_cell/PRDM9_used_peak_tss_2K.bed")
DF_used=DF[as.character(Finally_used_peak$V1),]

DSB <- c("PreL","L", "Z","PreP")
projHeme6 <- addTrajectory(
                                        ArchRProj = projHeme6,
                                        name = "DSB",
                                        groupBy = "Cell_type",
                                        preFilterQuantile = 0.99999,
                                        postFilterQuantile = 0.99999,
                                        trajectory = DSB,
                                        reducedDims = "IterativeLSI",
                                        #dof = 500,
                                        embedding = "UMAP",
                                        force = TRUE
                                                )

meta_L_Z=as.data.frame(projHeme6@cellColData[projHeme6@cellColData$Cell_type%in%c("L", "Z"),])
meta_L_Z$DSB_num=colSums(DF_used[,row.names(meta_L_Z)])
meta_L_Z=meta_L_Z[which(!is.na(meta_L_Z$DSB)),]
meta_L_Z=meta_L_Z[order(meta_L_Z$DSB),]
meta_L_Z$ID=seq(1:dim(meta_L_Z)[1])
meta_L_Z=na.omit(meta_L_Z)


######  DSB pattern barplot  ######

ggplot(data=meta_L_Z,mapping=aes(x=ID,y=DSB_num))+
        geom_bar(stat='identity',mapping=aes(fill=Clusters2))+
        #geom_smooth(method='gam',colour='blue',alpha=0.1,se = F)+
        scale_fill_manual(
                                        labels=c('L','Z'),
                                values=c('#7ec244','#c92c35'))+
    theme_classic()



######  DSB pattern heatmap  ######

data_df=DF_used[,row.names(meta_L_Z)]
Row_ann=meta_L_Z[,c("Cell_type","Sample")]
cols_23=c("#5cc2c4","#f5a9a6","#223a58","#f8e1a0","#f15e52","#007e70","#fa8428","#f91964","#dbc976",
                  "#f93faa","#074ba9","#f8d703","#8791a6","#c6156a","#4a211b","#442871","#989602","#ddfc00","#00ffd4",
                  "#0835cf","#ff0000","#00ff00","#0013ff")
ann.color <- list(Sample = setNames(ArchRPalettes$stallion[1:length(levels(Row_ann$Sample))],levels(Row_ann$Sample)),
                                                    Cell_type = setNames(c("#81bf2c","#c92327"),c("L","Z")),
                                                                        Chr=setNames(cols_23, as.character(unique(col_ann$Chr))  ))

 p=pheatmap(data_df,
                    show_rownames = F,show_colnames = F,
                        cluster_cols = FALSE,
                        clustering_method = "ward.D2",
                        annotation_colors = ann.color,
                        annotation_row = col_ann,
                        annotation_col = Row_ann,
                        cutree_rows=4,
                        color = paletteContinuous(set = "whiteRed", n = 100,reverse = FALSE))

ggsave("/date/xiehaoling/Sperm_atac/OA/PRDM9/heatmap_ALL_cell/L_Z_test_4_20_2K.png",p,width = 10, height = 13, dpi = 300)

row_cluster <- as.data.frame(cutree(p$tree_row,k=4))
write.csv(row_cluster,"L_Z_pheatmap_row_cluster_4_2K.csv")

######  Enrichment analysis of four types of peaks   ######
library("LOLA")
dbPath = system.file("extdata", "hg19", package="LOLA")
regionDB = loadRegionDB(dbPath)

peak_cluster1 <- readPeakFile("PRDM9_peak_cluster_A_2K.bed")
peak_cluster2 <- readPeakFile("PRDM9_peak_cluster_B_2K.bed")
peak_cluster3 <- readPeakFile("PRDM9_peak_cluster_C_2K.bed")
peak_cluster4 <- readPeakFile("PRDM9_peak_cluster_D_2K.bed")
userUniverse <- readPeakFile("All_bg_Markers_peak.bed")

locResults = runLOLA(peak_cluster1, userUniverse, regionDB, cores=1)  # peak_cluster2,peak_cluster3,peak_cluster4

