library(ArchR)
library('magrittr')
set.seed(1)

inputFiles=getInputFiles("Path to all OA and NOA fragment files")
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


saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-OA_NOA_cell", load = FALSE)





projHeme1 <- addGroupCoverages(ArchRProj =  projHeme1, groupBy = "Clusters",force = TRUE)
projHeme2 <- addReproduciblePeakSet(
                                ArchRProj =  projHeme1,
                                groupBy = "Clusters",
                                peaksPerCell = 900,
                                maxPeaks = 700000,
                                excludeChr = c("chrM"),
                                pathToMacs2 = "/home/xiehaoling/.pyenv/versions/anaconda3-5.3.1/envs/macs2/bin/macs2"
                                        )

OA_NOA_projHeme <- addPeakMatrix( projHeme2)
OA_NOA_projHeme <- addIterativeLSI(
                                ArchRProj = OA_NOA_projHeme,
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

OA_NOA_projHeme <- addClusters(
                                input = OA_NOA_projHeme,
                                force = TRUE,
                                reducedDims = "IterativeLSI",
                                method = "Seurat",
                                name = "Clusters",
                                maxClusters = 100,
                                resolution = 7 )


OA_NOA_projHeme <- addUMAP(
                                ArchRProj = OA_NOA_projHeme,
                                force = TRUE,
                                reducedDims = "IterativeLSI",
                                name = "UMAP",
                                nNeighbors = 40,
                                minDist = 0.5,
                                metric = "cosine"
                                        )

addGeneScoreMatrix(
                                input = OA_NOA_projHeme,
                                genes = getGenes(OA_NOA_projHeme),
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
                                blacklist = getBlacklist(OA_NOA_projHeme),
                                threads = getArchRThreads(),
                                parallelParam = NULL,
                                subThreading = TRUE,
                                force = TRUE,
                                logFile = createLogFile("addGeneScoreMatrix")
                                )

OA_NOA_projHeme <- addMotifAnnotations(ArchRProj = OA_NOA_projHeme, motifSet = "cisbp", name = "Motif")
OA_NOA_projHeme <- addBgdPeaks(OA_NOA_projHeme)

OA_NOA_projHeme <- addDeviationsMatrix(
                                ArchRProj = OA_NOA_projHeme,
                                peakAnnotation = "Motif",
                                force = TRUE
                                )

saveArchRProject(ArchRProj = OA_NOA_projHeme , load = FALSE)


######   Identification of cell types   #######

OA_meta=read.csv("/date/xiehaoling/Sperm_atac/OA_Y/OA_cell_sample_infor.csv",row.names = 1)

j=1
OA_lable=c()
for (i in as.character(row.names(OA_NOA_projHeme@cellColData)))
            {
                        if (i %in% row.names(OA_meta))
                        {
                            OA_lable[j]=as.character(OA_meta[i,"Cell_type"])
                j=j+1
                    }else
                        {
                            OA_lable[j]="NOA"
                        j=j+1
                        }
                }

OA_NOA_projHeme@cellColData$OA_lable=OA_lable
OA_NOA_projHeme@cellColData$new_lable=OA_lable

for (i in as.data.frame(table(OA_NOA_projHeme@cellColData$Clusters))$Var1 %>% as.character())
            {

                            tmp_table=as.data.frame(table(OA_NOA_projHeme@cellColData[OA_NOA_projHeme@cellColData$Clusters==i,"OA_lable"]))

    if("NOA"%in%tmp_table$Var1)
                        {
                                max_num=tmp_table[tmp_table$Var1!="NOA",]$Freq %>% which.max()
                    max_lable=tmp_table[tmp_table$Var1!="NOA",]$Var1[max_num] %>% as.character()
                            OA_NOA_projHeme@cellColData[OA_NOA_projHeme@cellColData$Clusters==i & OA_NOA_projHeme@cellColData$OA_lable=="NOA","new_lable"]=max_lable

                            }

            }

OA_NOA_projHeme@cellColData$Cell_type=OA_NOA_projHeme@cellColData$new_lable



######   Developmental arrest    ######

for (i %in% c("NOA1","NOA2","NOA3","NOA4","NOA5","NOA6","NOA7","NOA8","NOA9","NOA10","NOA11","NOA12","NOA13","NOA14","NOA15"))
{
pdf(file = paste("/date/xiehaoling/Sperm_atac/OA_NOA/Plots/NOA_sample/",i,"_umap.pdf",sep =""), height=8, width=8)

pm3 <- plotEmbedding(ArchRProj = OA_NOA_projHeme, colorBy = "cellColData",rastr = FALSE,size=0.5, name = "Cell_type",
                                                    mbedding = "UMAP",labelMeans=FALSE,
                                                        highlightCells = row.names(OA_NOA_projHeme@cellColData[OA_NOA_projHeme@cellColData$Sample==i,])
                                                                                                      )
pm3
dev.off()
}

######   Differential gene activation in patients with AZFc microdeletions  ######
for (Azf_patients %in% c("NOA6_Z","NOA12_Z","NOA14_Z"))
{

markersPeaks <- getMarkerFeatures(
                                                ArchRProj = OA_NOA_projHeme,
                                                useMatrix = "GeneScoreMatrix",
                                                groupBy = "pa_diff_gene",
                                                bias = c("TSSEnrichment", "log10(nFrags)"),
                                                testMethod = "wilcoxon",
                                                bgdGroups = "O_Z",
                                                useGroups = Azf_patients
                                                                                                       )
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.5 & abs(Log2FC) >= 0.5")
write.csv(markerList[[1]],paste("/date/xiehaoling/Sperm_atac/Save-OA_NOA_chrY_cell/Patient_Type/diff_gene/OA_vs_",Azf_patients,"diff_gene.csv",sep = "_"))

}

######  Enrichment analysis of differential peaks during the Z period (by patient type)  ######
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

OA_NOA_projHeme <- addPeakAnnotations(ArchRProj = OA_NOA_projHeme, regions = ELE, name = "ELE")


markerTest <- getMarkerFeatures(
                                                    ArchRProj = OA_NOA_projHeme,
                                                        useMatrix = "PeakMatrix",
                                                        groupBy = "Patient_Type_celltype",
                                                        testMethod = "wilcoxon",
                                                        bias = c("TSSEnrichment", "log10(nFrags)"),
                                                        useGroups = "Type1_Z", #Type2_Z,Type3_Z
                                                        bgdGroups = "OA_Z"
                                                                                  )

markerList <- getMarkers(markerTest, cutOff = "FDR < 0.05 & abs(Log2FC) >= 0.5")


enrichRepeat_up <- peakAnnoEnrichment(
                                                        seMarker = markerTest,
                                                        ArchRProj = OA_NOA_projHeme,
                                                        peakAnnotation = "ELE",
                                                        cutOff = "FDR <= 0.05 & Log2FC >= 1"
                                                                        )

######  Differential gene activation during the Z period (by patient type)  ######

markersGS <- getMarkerFeatures(
                                                        ArchRProj = OA_NOA_projHeme,
                                                        useMatrix = "GeneScoreMatrix",
                                                        groupBy = "Patient_Type_celltype",
                                                        bgdGroups =c("Type1_Z","Type2_Z","Type3_Z","OA_Z"),
                                                        useGroups =c("Type1_Z","Type2_Z","Type3_Z","OA_Z"),
                                                        bias = c("TSSEnrichment", "log10(nFrags)"),
                                                        testMethod = "wilcoxon"
                                                                                           )


######  Differential motif activation during the Z period (by patient type)   ############

diffMotif <- getMarkerFeatures(
                                                        ArchRProj = OA_NOA_projHeme,
                                                        testMethod = "wilcoxon",
                                                        useGroups = "OA_Z",
                                                        bgdGroups = "Type1_Z", #Type2_Z,Type3_Z
                                                        binarize = FALSE,
                                                        useMatrix = "MotifMatrix",
                                                        groupBy = "Patient_Type_celltype",
                                                        k = 500,
                                                        bufferRatio = 1,
                                                        useSeqnames="z")

