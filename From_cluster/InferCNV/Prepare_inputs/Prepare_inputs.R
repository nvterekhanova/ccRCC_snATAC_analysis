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

library(EnsDb.Hsapiens.v86)

sample="CPT0020120013"

atac=readRDS(paste('../../2.Cell_annotation/Annot_20201027/$sampleID\_ATAC_cellTyped.rds',
sep=""))

Idents(atac)=atac$predicted.id
orig=as.data.frame(Idents(atac))
colnames(orig)=c('cell_type')
orig$cell_type=gsub(" ", "_",orig$cell_type)
orig$cell_type=ifelse(orig$cell_type=="Tumor","Tumor","Obs")

write.table(orig, paste(sample,'.Barcode_Annotation.txt',sep=""), col.names=FALSE,sep="\t",quote=FALSE)

DefaultAssay(atac)='ACTIVITY'
#mat=atac@assays$RNA
mat=GetAssayData(atac, assay='ACTIVITY')
mat=as.data.frame(mat)
write.table(mat,paste(sample,'.RNA_Count.tsv',sep=''),sep='\t',quote=FALSE)
