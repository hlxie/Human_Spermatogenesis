library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
#library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(Rcpp)
library(viridis)

set.seed(1)


#----------------------------
#Reads count stats for a 5M window per cell
#Adapted form https://github.com/GreenleafLab/10x-scATAC-2019/blob/master/code/08_Run_scCNV_v2.R
#----------------------------

"%ni%" <- Negate("%in%")

countInsertions <- function(query, fragments, by = "RG"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1],
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)),
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

makeWindows <- function(genome, blacklist, windowSize = 5e6, slidingSize = 5e6){
        chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
        chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
        windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
        mcols(windows)$wSeq <- as.character(seqnames(windows))
        mcols(windows)$wStart <- start(windows)
        mcols(windows)$wEnd <- end(windows)
        message("Subtracting Blacklist...")
        windowsBL <- lapply(seq_along(windows), function(x){
                        if(x %% 100 == 0){
                                message(sprintf("%s of %s", x, length(windows)))
                        }
                        gr <- GenomicRanges::setdiff(windows[x,], blacklist)
                        mcols(gr) <- mcols(windows[x,])
                        return(gr)
                })
        names(windowsBL) <- paste0("w",seq_along(windowsBL))
        windowsBL <- unlist(GenomicRangesList(windowsBL), use.names = TRUE)
        mcols(windowsBL)$name <- names(windowsBL)
        message("Adding Nucleotide Information...")
        windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
        windowNuc <- lapply(seq_along(windowSplit), function(x){
                message(sprintf("%s of %s", x, length(windowSplit)))
            chrSeq <- Biostrings::getSeq(genome,chromSizes[which(seqnames(chromSizes)==names(windowSplit)[x])])
            grx <- windowSplit[[x]]
            aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
            mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
            mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
            return(grx)
          }) %>% GenomicRangesList %>% unlist %>% sortSeqlevels %>% sort
        windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
        windowNuc
}

scCNA <- function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrM")){

        #Keep only regions in filtered chromosomes
        windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
        fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
        windows <- windows[seqnames(windows) %ni% remove]
        fragments <- fragments[seqnames(fragments) %ni% remove]

        #Count Insertions in windows
        message("Getting Counts...")
        counts <- countInsertions(windows, fragments, by = "RG")[[1]]
        message("Summarizing...")
        windowSummary <- GenomicRangesList()
        countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
        for(x in seq_along(unique(mcols(windows)$name))){
                if(x %% 100 == 0){
                        message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
                }
                idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
                wx <- windows[idx,]
                wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
                mcols(wo)$name <- mcols(wx)$name[1]
                mcols(wo)$effectiveLength <- sum(width(wx))
                mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
                mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
                mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
                mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
                countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
                windowSummary[[x]] <- wo
        }
        windowSummary <- unlist(windowSummary)

        #Keep only regions with less than 0.1% N
        keep <- which(windowSummary$N < 0.001)
        windowSummary <- windowSummary[keep,]
        countSummary <- countSummary[keep,]

        #Now determine the nearest neighbors by GC content
        message("Computing Background...")
        bdgMean <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        bdgSd <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        log2FC <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        z <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
        pval <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))

        for(x in seq_len(nrow(countSummary))){
                if(x %% 100 == 0){
                        message(sprintf("%s of %s", x, nrow(countSummary)))
                }
                #Get Nearest Indices
                idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
                idxNN <- idxNN[idxNN %ni% x]
                #Background
                if(any(colMeans(countSummary[idxNN, ])==0)){
                        if(force){
                                message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
                        }else{
                                stop("Background Mean = 0!")
                        }
                }
                bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
                bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
                log2FC[x, ] <- log2((countSummary[x, ]+1e-5) / (bdgMean[x, ]+1e-5))
                z[x, ] <- (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
                pval[x, ] <- 2*pnorm(-abs(z[x, ]))
        }
        padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
        CNA <- matrix(0, nrow=nrow(countSummary), ncol=ncol(countSummary))
        CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

        se <- SummarizedExperiment(
                assays = SimpleList(
                                CNA = CNA,
                                counts = countSummary,
                                log2FC = log2FC,
                                padj = padj,
                                pval = pval,
                                z = z,
                                bdgMean = bdgMean,
                                bdgSd = bdgSd
                        ),
                rowRanges = windowSummary
        )
        colnames(se) <- colnames(counts)

        return(se)
}


#----------------------------
# Get inputs,using NOA10 as an example
#----------------------------
blacklist <- import.bed("hg19-blacklist.v2.bed")
windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg19, blacklist = blacklist)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
minFrags <- 500
filterFrags <- 1500
filterTSS <- 4
file_fragments <- "NOA10.fragments.tsv.gz"
fragments <- data.frame(readr::read_tsv(file_fragments, col_names=FALSE))
fragments <- GRanges(
                                seqnames = fragments[,1],
                                IRanges(fragments[,2]+1, fragments[,3]),
                                RG = fragments[,4],
                        N = fragments[,5]
                         )

tabRG <- table(fragments$RG)
keep <- names(tabRG)[which(tabRG >= minFrags)]
fragments <- fragments[fragments$RG %in% keep,]
fragments <- sort(sortSeqlevels(fragments))

cnaObj <- scCNA(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrM"))
saveRDS(cnaObj, "NOA10_GC.rds")

#----------------------------
#Calculate aneuploid chromosomes compared to the OA control
#----------------------------

library(dplyr)
library(DNAcopy)
library(ggplot2)
library(pheatmap)
library(gplots)
library(DNAcopy)

Nor_data=read.csv("Nor_data.csv") #load control file
meta=read.csv("/date/xiehaoling/Sperm_atac/OA_NOA_Y/OA_NOA_cell_sample_infor.csv",row.names = 1)

name=c()
j=1
for(i in row.names(meta))
{
        name[j]=strsplit(i,"#")[[1]][2]
        j=j+1

}

meta$name=name
meta_plot=data.frame(name=meta$name,Clusters2=meta$Clusters2)
row.names(meta_plot)=meta_plot$name



NOA10_GC=readRDS("NOA10_GC.rds")
cellname_NOA10=NOA10_GC@colData@rownames
#data_NOA10=NOA10_GC@assays@data[["counts"]]
data_NOA10=NOA10_GC@assays@data$counts

colnames(data_NOA10)=NOA10_GC@colData@rownames

window=as.data.frame(NOA10_GC@rowRanges)
window$num=window$end/10000000

sort_data=meta_plot[colnames(data_NOA10),]
sort_data$Clusters2=as.character(sort_data$Clusters2)
sort_data=na.omit(sort_data)

vec <- c("SSC_1", "SSC_2","Diff_ing_SPG", "Diff_ed_SPG", "PreL", "L","Z", "PreP", "P", "D", "Second","S1", "S2", "S3")

new_dataset <- left_join(data.frame(Clusters2 = vec),sort_data,by = "Clusters2")
new_dataset=na.omit(new_dataset)

nor_data_NOA10=data_NOA10
for(i in 1:dim(data_NOA10)[2])
        nor_data_NOA10[,i]=data_NOA10[,i]/colSums(data_NOA10)[i]
    NOA10_CNV=as.data.frame(matrix(nrow=489, ncol=nrow(new_dataset)))
        NOA10_CNV_log = as.data.frame(matrix(nrow=489, ncol=nrow(new_dataset)))
        for (i in 1:dim(new_dataset)[1] )
                   {
                         cell_name=as.character(new_dataset[i,2])
                 cell_type=as.character(new_dataset[i,1])
                 NOA10_CNV[,i]=nor_data_NOA10[,cell_name]/Nor_data[,cell_type]
                 NOA10_CNV_log[,i]=log2(nor_data_NOA10[,cell_name]/Nor_data[,cell_type])

                  }
colnames(NOA10_CNV)=as.character(new_dataset$name)
colnames(NOA10_CNV_log)=as.character(new_dataset$name)


NOA10_CNV_log[NOA10_CNV_log< -10]= -10

window_number=c()
j=1
for (i in window$seqnames)
        {
          window_number[j]=strsplit(i,"r")[[1]][2]
      j=j+1

         }

window_number=gsub("X",23,window_number)
window_number=gsub("Y",24,window_number)
window$window_number=as.numeric(window_number)


CNV_cell=list()
j=1
NOA10_CNV_Seg=NOA10_CNV_log
for (i in colnames(NOA10_CNV_log))
{
            test1 <- segment(CNA(NOA10_CNV_log[,i], sampleid=i,window$window_number, window$num),min.width=2)
                CNV_cell[[j]]=na.omit(test1$output[abs(test1$output$seg.mean)>2,])
                j=j+1

                n=1
                llist=list()
                test1$output$num.mark=test1$segRows$endRow-test1$segRows$startRow+1

                for (m in 1:dim(test1$output)[1])
                        {
                          llist[[n]]=rep(test1$output[m,"seg.mean"],test1$output[m,"num.mark"])
                          n=n+1

                         }
                NOA10_CNV_Seg[,i]=unlist(llist)
}
NOA10_CNV_ALL=Reduce(function(x, y) merge(x, y, all= TRUE ), CNV_cell)
cell=NOA10_CNV_ALL[abs(NOA10_CNV_ALL$seg.mean)>3,]$ID
new_dataset=unique(new_dataset)

row.names(new_dataset)=as.character(new_dataset$name)
NOA10_CNV_ALL$cell_type=new_dataset[gsub("[.]","-",NOA10_CNV_ALL$ID),"Clusters2"]

write.csv(NOA10_CNV_ALL,"NOA10_CNV_ALL.csv")
write.csv(NOA10_CNV,"NOA10_CNV.csv")
write.csv(NOA10_CNV_Seg,"NOA10_CNV_Seg.csv")
write.csv(NOA10_CNV_log,"NOA10_CNV_log.csv")


#----------------------
#Statistical analysis of sperm with aneuploidy
#----------------------

CNV_ALL=read.csv("NOA10_CNV_ALL.csv")
CNV_ALL$X=NULL
CNV=read.csv("NOA10_CNV.csv")
CNV$X=NULL
CNV_Seg=read.csv("NOA10_CNV_Seg.csv")
CNV_Seg$X=NULL
CNV_log=read.csv("NOA10_CNV_log.csv")
CNV_log$X=NULL

meta=read.csv("/date/xiehaoling/Sperm_atac/OA_NOA_Y/OA_NOA_cell_sample_infor.csv",row.names = 1)
name=c()
j=1
for(i in row.names(meta))
        {
          name[j]=strsplit(i,"#")[[1]][2]
          j=j+1
        }

meta$name=name
meta$name2=gsub("-",".",meta$name)
row.names(meta)=meta$name2
total_cell=dim( meta[meta$Sample=="NOA10" & meta$Clusters2%in%c("Second","S1","S2","S3"),] )[1]

window=read.csv("/date/xiehaoling/Sperm_atac/CNV/window.csv")

CNV[CNV<0]=0
CNV[CNV>5]=5

CNV$chr=window$seqnames
CNV$sub_chr=window$seqnames_sub

com_cell_num= dim(CNV)[2]-2

group_mean = matrix(ncol = com_cell_num, nrow = 34)
for (i in 1:(dim(CNV)[2]-2))
        {
          group_mean[,i] <- aggregate(CNV[,i], list(CNV$sub_chr), mean)$x


        }
colnames(group_mean)=colnames(CNV)[1:com_cell_num]
row.names(group_mean)=aggregate(CNV[,i], list(CNV$sub_chr), mean)$Group.1

######CNV_do

chr=row.names(which(group_mean < 0.2, arr.ind=TRUE))
CNV_do=as.data.frame(which(group_mean < 0.2, arr.ind=TRUE))
CNV_do$chr=chr
CNV_do$name=colnames(group_mean[,CNV_do$col])
CNV_do$type=meta[CNV_do$name,"Clusters2"]
CNV_do$CNV=rep("do",dim(CNV_do)[1])
CNV_do=na.omit(CNV_do)

autosomal_do=CNV_do[CNV_do$row>=1 & CNV_do$row<=31 & CNV_do$chr!="un" & CNV_do$type%in%c("Second","S1","S2","S3"),]

CNV_do_XY=CNV_do[CNV_do$type%in%c("S1","S2","S3","S4") & CNV_do$chr%in%c("chrX","chrY"),]
cell1=as.character(as.data.frame(table(CNV_do_XY$name))[ as.data.frame(table(CNV_do_XY$name))$Freq==2,"Var1" ])
XY_do=CNV_do[CNV_do$name%in%c(cell1),]

tmp1=rbind(autosomal_do,XY_do)
write.csv(rbind(tmp1,c(total_cell,dim(autosomal_do)[1], length(cell1),"total_cell","autosomal_do","XY_do" )),"./NOA10_CNV_cell_do")



#######CNV_up
chr=row.names(which(group_mean >3, arr.ind=TRUE))
CNV_up=as.data.frame(which(group_mean >3, arr.ind=TRUE))

CNV_up$chr=chr
CNV_up$name=colnames(group_mean[,CNV_up$col])
CNV_up$type=meta[CNV_up$name,"Clusters2"]
CNV_up$CNV=rep("up",dim(CNV_up)[1])

CNV_up=na.omit(CNV_up)
up=CNV_up[ CNV_up$chr!="un" & CNV_up$type%in%c("Second","S1","S2","S3"),]

 ######CNV_XY
chr=row.names(which(group_mean >=0.5, arr.ind=TRUE))
CNV_XY=as.data.frame(which(group_mean >=0.5, arr.ind=TRUE))
CNV_XY$chr=chr
CNV_XY$name=colnames(group_mean[,CNV_XY$col])
CNV_XY$type=meta[CNV_XY$name,"Clusters2"]
CNV_XY$CNV=rep("up",dim(CNV_XY)[1])

ss=as.data.frame(table(CNV_XY[CNV_XY$type%in%c("Second","S1","S2","S3") & CNV_XY$row%in%c(33,32),"name"]))
cell2=as.character(ss[ss$Freq>1,"Var1"])

CNV_XY_data=CNV_XY[CNV_XY$name%in%c(cell2) & CNV_XY$row%in%c(33,32) ,]
tmp1=rbind(up,CNV_XY_data)

write.csv(rbind(tmp1,c(total_cell,dim(up)[1],length(cell2),"total_cell","up","XY_up")),"./NOA10_CNV_cell_up")
                    
                     
                      
