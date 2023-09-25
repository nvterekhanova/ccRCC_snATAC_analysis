###Need to map clusters:
#system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

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

library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers =5)
options(future.globals.maxSize = 10000 * 1024^2)


ATAC=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'cicero.v3.20210725.rds',sep=''))

cell_t=read.table('../../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type

genes=c("ABCC3", "ABLIM3",  "CA9","COL23A1", "CP", "EGFR", "ENPP3", "EPHA6", "FTO","KCTD3","NDRG1", 
 "PCSK6","PHKA2", "PLEKHA1", "PLIN2", "SEMA6A","SHISA9", "SLC6A3", "SNAP25", "TGFA", "UBE2D2")

####First, calculate FChange for Tumor vs PT:
ATAC$test=ATAC$cell_type
ATAC$test=ifelse(ATAC$cell_type %in% c('Tumor','EMT tumor cells'),'Tumor',ATAC$test)
Idents(ATAC)=ATAC$test
DefaultAssay(ATAC)='ATACGeneActivity'
ident.use.1='Tumor'
ident.use.2='PT'
fc.results <- FoldChange(
  object = ATAC,
  assay="ATACGeneActivity",
  slot = "data",
  ident.1 = ident.use.1,
  ident.2 = ident.use.2,
  base = 2
)

fc.results$Gene=row.names(fc.results)
write.table(fc.results,'out/FoldChange_ATACGeneActivity_Tumor_vs_PT.20210924.tsv',sep='\t',
row.names=F,quote=F)


###Second, calculate FChange for Tumor vs All Other cells:

ATAC$test=ATAC$cell_type
ATAC$test=ifelse(ATAC$cell_type %in% c('Tumor','EMT tumor cells'),'Tumor','Other')
Idents(ATAC)=ATAC$test
DefaultAssay(ATAC)='ATACGeneActivity'
ident.use.1='Tumor'
ident.use.2='Other'
fc.results <- FoldChange(
  object = ATAC,
  assay="ATACGeneActivity",
  slot = "data",
  ident.1 = ident.use.1,
  ident.2 = ident.use.2,
  base = 2
)

fc.results$Gene=row.names(fc.results)
write.table(fc.results,'out/FoldChange_ATACGeneActivity_Tumor_vs_AllOtherCells.20210924.tsv',sep='\t',
row.names=F,quote=F)



