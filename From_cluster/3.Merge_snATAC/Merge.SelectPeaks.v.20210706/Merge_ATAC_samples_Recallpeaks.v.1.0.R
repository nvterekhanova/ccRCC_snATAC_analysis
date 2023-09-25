library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM
date='20210706'

samples=c('')


tab=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/',
'snRNA_processed_Yige/meta_data.20200505.v1.tsv',sep=''),sep='\t',header=TRUE)
tab=tab[tab$Aliquot.snRNA %in% samples,]
samples=tab$Aliquot.snRNA
piece_ids=tab$Aliquot.snRNA.WU

atac=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    atac[[i]]=readRDS(paste("../../1.Create_rds/out/",samples[i],"/",samples[i],"_processed_atac.rds",
sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['RNA']]<-NULL
    atac[[i]][['peaks']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}

#####To obtain the best results - use ALL peaks!

combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#peaks.use=combined.peaks
####Now using MACS2-peak calling:
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#For testing purposes only:
#peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)

#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells are pass those filters in the tutorial.

matrix.counts=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
    ) 
}


for (i in 1:length(samples)){
atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],
fragments=Fragments(atac[[i]]@assays$X500peaksMACS2))
atac[[i]]$dataset=samples[i]
DefaultAssay(atac[[i]])<-'peaksinters'
###remove other assay
#atac[[i]][['X500peaksMACS2']]<-NULL
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
saveRDS(combined, paste('28_ccRCC_snATAC.',date,'.rds',sep=''))





