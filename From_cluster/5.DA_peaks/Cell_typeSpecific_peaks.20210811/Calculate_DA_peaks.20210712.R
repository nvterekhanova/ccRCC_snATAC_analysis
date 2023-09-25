system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

####FOR PBRM1-analysis -- INCLUDE 1287 sample in the non-PBRM1-group

	       
###For some samples (with many cells >6K) python-package used in chromVar doesn't work properly; Need to use this 2 commands:
###export OMP_NUM_THREADS=1
###export USE_SIMPLE_THREADED_LEVEL3=1
###export OPENBLAS_NUM_THREADS=1

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
#options(future.globals.maxSize = 100 * 1024^3) # for 200 Gb RAM
options(future.globals.maxSize = 10000 * 1024^2)


ATAC=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'20210707.rds',sep=''))


###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type

ATAC$test=ATAC$cell_type

ATAC$test=as.character(ATAC$cell_type)
ATAC$test=ifelse(ATAC$cell_type %in% c('DC','CD4+ T-cells','CD8+ T-cells','Macrophages','NK cells'),"Immune",
ATAC$test)
ATAC$test=ifelse(ATAC$cell_type %in% c('Endothelial cells','Fibroblasts'),"Stroma",ATAC$test)



peak.data <- GetAssayData(object = ATAC, assay = 'peaksMACS2', slot = "counts")
total_fragments_cell <- ATAC$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
ATAC <- AddMetaData(object = ATAC, metadata = frip, col.name = 'frip_500MACS2')
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_500MACS2')


tumor_piece_ids=unique(ATAC$Piece_ID)
cell_types=paste(tumor_piece_ids,"Tumor",sep='_')
###Check that cell type distribution looks OK:
table(ATAC@meta.data %>% select ('Piece_ID','cell_type'))

#Find Markers:

DefaultAssay(ATAC) <- 'peaksMACS2'
Idents(ATAC)=ATAC$test

#Now do the same using the DEGs-cutoffs for consistency  (from Yige):
#logfc.threshold.run <- 0
#min.pct.run <- 0.1
#min.diff.pct.run <- 0.1 ### change it to 0 similar to BAP1/PBRM1-peaks:


cell_types=c('Immune','Stroma','Tumor')
all_da_peaks=NULL
for (cell_type in cell_types){
da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_type,
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'peak_RF_500MACS2'
)
da_peaks$peak=rownames(da_peaks)
da_peaks$Sample=cell_type
all_da_peaks=rbind(all_da_peaks,da_peaks)
print(cell_type)
}
write.table(all_da_peaks, paste("out/da_peaks_cellGroup_Specific.minPct0.1.20210812.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)







