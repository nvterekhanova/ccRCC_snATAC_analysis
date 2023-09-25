#system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

###USE Peaaks that are more open/close in >=50% of BAP1-samples


library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyverse)

####First DOWN-peaks:
peaks=read_delim('data/EMT_vs_selectedEpithelialClusters.DEG_overlap_DAP.Consistent.20210928.v1.tsv',
delim='\t')
peaks=as.data.frame(peaks)

ATAC=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'cicero.v3.20210725.rds',sep=''))

peak.data <- GetAssayData(object = ATAC, assay = 'peaksMACS2', slot = "counts")
total_fragments_cell <- ATAC$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
ATAC <- AddMetaData(object = ATAC, metadata = frip, col.name = 'frip_500MACS2')
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_500MACS2')

barcode_map=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/',
'12.Transfer_tumorClusters_fromRNA/TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/',
'RNA_clusters_Tumor_PT_clusters.snATAC.20210917.tsv',sep=''),sep='\t',header=T)
rownames(barcode_map)=barcode_map$sample_barcode

atac=ATAC
ATAC=subset(atac, cells=barcode_map$sample_barcode)

orig_1=as.data.frame(ATAC$orig.ident)
barcode_map=barcode_map[rownames(orig_1),]
ATAC$cell_group_ID=barcode_map$cell_group
Idents(ATAC)=ATAC$cell_group_ID


x=ATAC@assays$peaksMACS2
x1=ATAC@assays$peaksMACS2
x1=x1[rownames(x1) %in% peaks$peak,]
x1=as.data.frame(x1)
y=t(x1)
y=as.data.frame(y)

y$cell_group_ID=ATAC$cell_group_ID
#y$Case=ATAC$Piece_ID
#y$Cell_type=ATAC$cell_type

m=colSums(x)
y1=y
y1[,1:(ncol(y1)-1)]=y1[,1:(ncol(y1)-1)]/m
y1[,1:(ncol(y1)-1)]=y1[,1:(ncol(y1)-1)]*100000


###Now aggregate:
t1=aggregate(y1[,1:(ncol(y1)-1)],by=list(y1$cell_group_ID),FUN='mean',na.action = na.omit)

write.table(t1,"out/Accessibility.EMT_vs_selectedEpithelialClusters.DEG_overlap_DAP.Consistent.20210928.tsv",
sep='\t',quote=F,row.names=F)
