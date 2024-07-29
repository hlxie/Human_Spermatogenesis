library(ArchR)
set.seed(1)

inputFiles=getInputFiles("Path to all OA fragment files")
addArchRGenome("hg19")
addArchRThreads(threads = 10)

ArrowFiles <- createArrowFiles(
                inputFiles = inputFiles,
                sampleNames = names(inputFiles),
                filterTSS = 4, #Dont set this too high because you can always increase later
                                excludeChr = c("chrM"),
                filterFrags = 10000,
                maxFrags = 1000000000,
                addTileMat = TRUE,
                addGeneScoreMat = FALSE
                )


projHeme1 <- ArchRProject(
                            ArrowFiles = ArrowFiles,
                            outputDirectory = "OA_NOA_cell",
                            copyArrows = FALSE
                            )


projHeme1 <- addIterativeLSI(
                ArchRProj = projHeme1,
                force = TRUE,
                sampleCellsPre = 15000,
                useMatrix = "TileMatrix",
                name = "IterativeLSI",
                iterations = 4,

                UMAPParams = list(
                n_neighbors = 40,
                min_dist = 0.6,
                metric = "cosine",
                verbose = FALSE,
                fast_sgd = TRUE
                  ),
                clusterParams = list(
                resolution = c(0.2),
                sampleCells = 30000,
                n.start = 10
                  ),
                totalFeatures = 500000,
                filterQuantile = 0.995,
                varFeatures = 35000,
                dimsToUse = 1:30,
                scaleTo = 100000,
                nPlot = 30000
                 )

projHeme1 <- addClusters(
                input = projHeme1,
                force = TRUE,
                reducedDims = "IterativeLSI",
                method = "Seurat",
                name = "Clusters",
                resolution = 2
                )

projHeme1 <- addUMAP(
                ArchRProj = projHeme1,
                force = TRUE,
                reducedDims = "IterativeLSI",
                name = "UMAP",
                nNeighbors = 40,
                minDist = 0.5,
                metric = "cosine"
                )


saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-OA_cell", load = FALSE)





projHeme1 <- addGroupCoverages(ArchRProj =  projHeme1, groupBy = "Clusters",force = TRUE)
projHeme2 <- addReproduciblePeakSet(
                                ArchRProj =  projHeme1,
                                groupBy = "Clusters",
                                peaksPerCell = 900,
                                maxPeaks = 700000,
                                excludeChr = c("chrM"),
                                pathToMacs2 = "/home/xiehaoling/.pyenv/versions/anaconda3-5.3.1/envs/macs2/bin/macs2"
                                        )

OA_projHeme <- addPeakMatrix( projHeme2)
OA_projHeme <- addIterativeLSI(
                                ArchRProj = OA_projHeme,
                                force = TRUE,
                                sampleCellsPre = 10000,#迭代使用细胞数
                                useMatrix = "PeakMatrix",
                                name = "IterativeLSI",
                                iterations = 3,

                                UMAPParams = list(
                                n_neighbors = 40,
                                min_dist = 0.6,
                                metric = "cosine",
                                verbose = FALSE,
                                fast_sgd = TRUE
                                        ),

                                clusterParams = list(
                                resolution = c(0.2),
                                sampleCells = 30000,
                                n.start = 10
                                        ),

                                totalFeatures = 650000,
                                filterQuantile = 0.995,
                                varFeatures = 25000,
                                dimsToUse = 1:30,
                                scaleTo = 1000000,
                                nPlot = 30000
                                        )

OA_projHeme <- addClusters(
                                input = OA_projHeme,
                                force = TRUE,
                                reducedDims = "IterativeLSI",
                                method = "Seurat",
                                name = "Clusters",
                                resolution = 2 )


OA_projHeme <- addUMAP(
                                ArchRProj = OA_projHeme,
                                force = TRUE,
                                reducedDims = "IterativeLSI",
                                name = "UMAP",
                                nNeighbors = 40,
                                minDist = 0.5,
                                metric = "cosine"
                                        )


addGeneScoreMatrix(
                                input = OA_projHeme,
                                genes = getGenes(OA_projHeme),
                                geneModel = "exp(-abs(x)/5000) + exp(-1)",
                                matrixName = "GeneScoreMatrix",
                                extendUpstream = c(1000, 1e+05),
                                extendDownstream = c(1000, 1e+05),
                                geneUpstream = 5000,
                                geneDownstream = 0,
                                useGeneBoundaries = TRUE,
                                useTSS = FALSE,
                                extendTSS = FALSE,
                                tileSize = 500,
                                ceiling = 4,
                                geneScaleFactor = 5,
                                scaleTo = 10000,
                                excludeChr = c("chrM"),
                                blacklist = getBlacklist(OA_projHeme),
                                threads = getArchRThreads(),
                                parallelParam = NULL,
                                subThreading = TRUE,
                                force = TRUE,
                                logFile = createLogFile("addGeneScoreMatrix")
                                )

OA_projHeme <- addMotifAnnotations(ArchRProj = OA_projHeme, motifSet = "cisbp", name = "Motif")
OA_projHeme <- addBgdPeaks(OA_projHeme)

OA_projHeme <- addDeviationsMatrix(
                                ArchRProj = OA_projHeme,
                                peakAnnotation = "Motif",
                                force = TRUE
                                )

plotVarDev <- getVarDeviations(OA_projHeme, name = "MotifMatrix", plot = TRUE)

saveArchRProject(ArchRProj = OA_projHeme , load = FALSE)


############    QC plot    ################

library(ggpointdensity)

Pre_QC=readRDS("/date/xiehaoling/Sperm_atac/NOA/NOA_ArchR_10000/QC_ALL/OA_NOA-Pre-Filter-Metadata.rds")
Pre_QC_meta=as.data.frame(Pre_QC)

pdf(file = "/date/xiehaoling/Sperm_atac/NOA/NOA_ArchR_10000/QC_ALL/QC_NOA_scatter.pdf", height=8, width=8)
scatter <-ggplot(Pre_QC_meta,aes(log10(nFrags),TSSEnrichment))+
                     geom_pointdensity(adjust = 4)+
                     theme_bw()+
                             scale_color_distiller(palette = "Reds", direction = 1)+
                             geom_vline(xintercept=4,linetype='dashed', color='#8B1A1A', size=1)+
                                 geom_hline(yintercept=4,linetype='dashed', color='#8B1A1A', size=1)
scatter

dev.off()


###########    Identification of cell types    ########

markers <- c("UTF1", "REC8", "DPH7", "DMC1", "SYCP1", "SYCP3", "SYCP2", "PIWIL1","OVOL1","TXNDC2","PRM1","DMRT1",
                         "OVOL2","INSL3","DCN","CD3D","BEX1","CD14","PECAM1")

OA_projHeme <- addImputeWeights(OA_projHeme, reducedDims = "IterativeLSI")
p <- plotEmbedding(
                        ArchRProj = OA_projHeme,
                        colorBy = "GeneScoreMatrix",
                        name = markers,
                        embedding = "UMAP",
                        pal=paletteContinuous(set = "blueYellow", n = 100,reverse = FALSE),
                        imputeWeights = getImputeWeights(OA_projHeme),
                        size = 0.2,
                        quantCut = c(0.8, 0.95),
                        plotAs = "points",
                        rastr = TRUE
                        )

markersGS <- getMarkerFeatures(
                        ArchRProj = OA_projHeme,
                        useMatrix = "GeneScoreMatrix",
                        groupBy = "Cell_type",
                        bias = c("TSSEnrichment", "log10(nFrags)"),
                        testMethod = "wilcoxon"
                        )

heatmapGS <- markerHeatmap(
                        seMarker = markersGS,
                        cutOff = "FDR <= 0.05 & Log2FC >= 1",
                        labelMarkers = markers,
                        transpose = TRUE
                                )

###########   Identification of cell type-specific peaks   ############

markersPeaks <- getMarkerFeatures(
                        ArchRProj = OA_projHeme,
                        useMatrix = "PeakMatrix",
                        groupBy = "Cell_type",
                        useGroups=c("SSC_1","SSC_2","Diff_ing_SPG","Diff_ed_SPG","PreL","L","Z","PreP","P","Second","S1","S2","S3"),
                        bgdGroups=c("SSC_1","SSC_2","Diff_ing_SPG","Diff_ed_SPG","PreL","L","Z","PreP","P","Second","S1","S2","S3"),
                        bias = c("TSSEnrichment", "log10(nFrags)"),
                        testMethod = "wilcoxon"
                                )


heatmapPeaks <- markerHeatmap(
                        seMarker = markersPeaks,
                        cutOff = "FDR <= 0.05 & Log2FC >= 1",
                        transpose = TRUE
                                )


##########    cell type-specific peak enrichment analysis  #########

ELE <- c(
                 CpG = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.CGI.bed",
         Exon = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.Exon.bed",
                 Intergenic = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.Intergenic.bed",
                 Intragenic = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.Intragenic.bed",
                 Intron = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.Intron.bed",
                 Promoter = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.Promoter.bed",
                 Genebody = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.Genebody.bed",
                 LINE1 = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.L1.bed",
                 LINE2 = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.L2.bed",
                 LINE = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.LINE.bed",
                 LTR = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.LTR.bed",
                 SINE = "/datf/xiehaoling/reference/Database_Meth/hg19/sub_groups/XHL/hg19.SINE.bed"
         )

OA_projHeme <- addPeakAnnotations(ArchRProj = OA_projHeme, regions = ELE, name = "ELE")

markersPeaks_growth <- getMarkerFeatures(
                                                        ArchRProj = projHeme6,
                                                        useMatrix = "PeakMatrix",
                                                        groupBy = "Cell_type",
                                                        useGroups=c("SSC_1"),
                                                        bgdGroups=c("SSC_2","Diff_ing_SPG","Diff_ed_SPG","PreL","L","Z","PreP","P","Second","S1","S2","S3"),
                                                        bias = c("TSSEnrichment", "log10(nFrags)"),
                                                        testMethod = "wilcoxon"
                                                                        )

enrichRepeat_up <- peakAnnoEnrichment(
                                                        seMarker = markersPeaks_growth,
                                                        ArchRProj = projHeme6,
                                                        peakAnnotation = "ELE",
                                                        cutOff = "FDR <= 0.05 & Log2FC >= 1"
                                                                        )


########   Trajectory  analysis   #########

Sperm=c("SSC_1","SSC_1","Diff_ing_SPG","Diff_ed_SPG","PreL","L","Z","PreP","P","D","Second","S1","S2","S3")

OA_projHeme <- addTrajectory(
                        ArchRProj = OA_projHeme,
                        name = "Sperm",
                        groupBy = "Cell_type",
                        preFilterQuantile = 0.99,
                        postFilterQuantile = 0.99,
                        trajectory = Sperm,
                        reducedDims = "IterativeLSI",
                        #dof = 500,
                        embedding = "UMAP",
                        force = TRUE
                        )


gene=c('RFX2','USF2','CTCF','FOS','JUN','JUND','YY1','MEF2B','EBF1','KLF14','MZF1','SP6','SP9','KLF16','KLF4','SP1')
motif=c('z:RFX2_724','z:USF2_26','z:CTCF_177','z:FOS_137','z:JUN_143','z:JUND_124','z:YY1_173','z:MEF2B_643','z:EBF1_67','z:KLF14_251','z:MZF1_171','z:SP6_275','z:SP9_283','z:KLF16_205','z:KLF4_208','z:SP1_267')

pdf(file = "Marker_Trajectory.pdf", height=4, width=16)

for(i in 1:16)
            {
        p2 <- plotTrajectory(OA_projHeme, trajectory = "Sperm", colorBy = "MotifMatrix", name = motif[i], size = 0.1, continuousSet = "solarExtra")
    p1 <- plotTrajectory(OA_projHeme, trajectory = "Sperm", colorBy = "GeneScoreMatrix", name = gene[i], size = 0.1, continuousSet = "blueYellow")
        print(  ggAlignPlots(p1[[2]], p2[[2]], type = "h")  )

            }

dev.off()


trajmotif <- getTrajectory(ArchRProj = OA_projHeme, name = "Sperm", useMatrix = "MotifMatrix", log2Norm = F)
trajPeak  <- getTrajectory(ArchRProj = OA_projHeme, name = "Sperm", useMatrix = "PeakMatrix", log2Norm = F)
trajGSM  <- getTrajectory(ArchRProj = OA_projHeme, name = "Sperm", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)

corGSM_MM <- correlateTrajectories(trajGSM, trajmotif,corCutOff = 0.4, varCutOff1 = 0.5, varCutOff2 = 0.5,)

trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajmotif2 <- trajmotif[corGSM_MM[[1]]$name2, ]

trajCombined <- trajmotif2
assay(trajCombined) <- t(apply(assay(trajmotif2), 1, scale)) + t(apply(assay(trajGSM2 ), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

tTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "blueYellow"), varCutOff = 0.5,labelTop = 1,
                                                                rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajmotif2, pal = paletteContinuous(set = "comet"), varCutOff = 0.5,labelTop = 1,
                                                        labelMarkers =paste(rep("z:",dim(TF)[1]),as.character(TF$x),sep=""), rowOrder = rowOrder)

pdf(file = "TrajectoryHeatmap.pdf", height=7, width=8)

ht1 +ht2

dev.off()


#######   Identification of X/Y sperm  #######

library(AUCell)

GeneScoreMatrix=getMatrixFromProject(OA_projHeme,useMatrix = "GeneScoreMatrix")
GeneScoreMatrix_data=as.matrix(GeneScoreMatrix@assays@data$GeneScoreMatrix)
row.names(GeneScoreMatrix_data)=GeneScoreMatrix@elementMetadata$name

sperm_cell=row.names(OA_projHeme@cellColData[OA_projHeme@cellColData$Cell_type%in%c("Second","S1","S2","S3"),])

cells_rankings = AUCell_buildRankings(GeneScoreMatrix_data[,sperm_cell])

geneSets=list()
geneSets_chrX=toupper(as.character(chrX_Nomir_genename)) # chrX_Nomir_genename refers to genes on chrX, excluding MIR genes
geneSets[["chrX"]]=geneSets_chrX

cells_AUC=AUCell_calcAUC(geneSets,cells_rankings,aucMaxRank=nrow(cells_rankings)*0.5)
cells_ass= AUCell_exploreThresholds(cells_AUC,plotHist = TRUE,nCores = 1,assign=TRUE)

AUC_XY=as.data.frame(t(as.data.frame(getAUC(cells_AUC))))
AUC_XY$type=AUC_XY$chrX
AUC_XY$type[AUC_XY$type>0.09]="X"
AUC_XY$type[AUC_XY$type<0.07]="Y"
AUC_XY$type[(AUC_XY$type>=0.07) & (AUC_XY$type<=0.09)]="un"

