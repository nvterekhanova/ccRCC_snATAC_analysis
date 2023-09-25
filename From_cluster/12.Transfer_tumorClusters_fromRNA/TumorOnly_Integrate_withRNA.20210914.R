library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

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
library(RColorBrewer)

library(EnsDb.Hsapiens.v86)
args = commandArgs(trailingOnly=TRUE)

samples=c('')

for (sample in samples){
print(sample)
rna_paths=read.table('data_snRNA/Paths_TumorCellOnlyReclustered_SeuratObject.20210805.v1.tsv',
sep='\t',header=T)
rna_paths=rna_paths[rna_paths$Aliquot.snRNA %in% samples,]
path_rna=rna_paths$Path_katmai[rna_paths$Aliquot.snRNA==sample]

RNA=readRDS(path_rna)
ATAC=readRDS(paste("TumorOnly.RDS/",sample,"_TumorOnly.20210913.rds",sep=""))

case=rna_paths$Aliquot.snRNA.WU[rna_paths$Aliquot.snRNA==sample]
DefaultAssay(ATAC)='X500peaksMACS2'

#####Labelling cell-types for the snRNA

meta=read.table('data_snRNA/Barcode_TumorPTLOH_ByCluster.20210903.v1.tsv',sep='\t',header=T)
meta$Sample=gsub('(.*)_.*','\\1',meta$sample_barcode)
meta_sample=meta[meta$Sample==sample,]
meta_sample$Cluster_group=gsub('(.*)_.*','\\1',meta_sample$cell_group)
meta_sample=meta_sample[meta_sample$Cluster_group==case,]
meta_sample$Cluster_ID=gsub('(.*)_(.*)','\\2',meta_sample$cell_group)
rownames(meta_sample)=gsub('(.*)_(.*)','\\2',meta_sample$sample_barcode)

orig_1=as.data.frame(RNA$orig.ident)
meta_sample=meta_sample[rownames(orig_1),]

RNA$Cluster_ID=meta_sample$Cluster_ID


###Assigning cell-types:

DefaultAssay(ATAC) <- 'RNA'

ATAC <- NormalizeData(
  object = ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC$nCount_RNA)
)

transfer.anchors <- FindTransferAnchors(
  reference = RNA,
  query = ATAC,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA$Cluster_ID,
  weight.reduction = ATAC[['lsi']],
  dims=2:30
)

ATAC <- AddMetaData(object = ATAC, metadata = predicted.labels)
ATAC$Cluster_ID=ATAC$predicted.id
ATAC$Piece_ID=case

###Making plots:

plot1 <- DimPlot(
  object = RNA,
  group.by = 'Cluster_ID',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(case,' snRNA-seq',sep=""))

plot2 <- DimPlot(
  object = ATAC,
  group.by = 'Cluster_ID',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(case,' snATAC-seq',sep=""))

pdf(paste("TumorOnly_RNAClusters_mapped/plots/",sample,"_",case,"_snRNA_ClusterIDs_mapped_to_snATAC.pdf",sep=""),
height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()

atac_meta=ATAC@meta.data
write.table(atac_meta,paste("TumorOnly_RNAClusters_mapped/metadata/",sample,"_TumorOnly.snATAC.meta.data",
sep=""),sep='\t',quote=FALSE)
}

