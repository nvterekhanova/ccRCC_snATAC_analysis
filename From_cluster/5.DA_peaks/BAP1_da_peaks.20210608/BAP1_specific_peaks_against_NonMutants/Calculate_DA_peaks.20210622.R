system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

####only BAP1-mutants vs non-mutants
	       
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

ATAC=readRDS(paste('../../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'20210707.rds',sep=''))


###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../../../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type

ATAC$test=ATAC$cell_type
ATAC$test=ifelse(ATAC$test=="Tumor",paste(as.character(ATAC$Piece_ID),"Tumor",sep="_"),
ATAC$test)

#we decided to include both-mutants
meta=read.table('data/Sample_categories.20210503.txt',sep='\t',header=T)
meta=meta[meta$Aliquot.WU %in% unique(ATAC$Piece_ID),]
#bap1_s=meta$Aliquot.WU[meta$Category %in% c('BAP1-mutant')]
bap1_s=meta$Aliquot.WU[meta$Category %in% c('BAP1-mutant','Both mutated')]
bap1_s=bap1_s[bap1_s!='C3L-01287-T1']
non_bap1_s=meta$Aliquot.WU[meta$Category %in% c('Non-mutant')]




ATAC$test='Other'
ATAC$test=ifelse(ATAC$cell_type=="Tumor" & ATAC$Piece_ID %in% bap1_s,ATAC$Piece_ID,ATAC$test)
ATAC$test=ifelse(ATAC$cell_type=="Tumor" & ATAC$Piece_ID %in% non_bap1_s,"NOT_BAP1_mutant",ATAC$test)
#ATAC$test='Other'
#ATAC$test=ifelse(ATAC$cell_type=="Tumor" & ATAC$Piece_ID %in% pbrm1_s,ATAC$Piece_ID,ATAC$test)
#ATAC$test=ifelse(ATAC$cell_type=="Tumor" & ATAC$Piece_ID %in% non_pbrm1_s,"NOT_PBRM1_mutant",ATAC$test)


peak.data <- GetAssayData(object = ATAC, assay = 'peaksMACS2', slot = "counts")
total_fragments_cell <- ATAC$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
ATAC <- AddMetaData(object = ATAC, metadata = frip, col.name = 'frip_500MACS2')
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_500MACS2')


tumor_piece_ids=unique(ATAC$Piece_ID)
tumor_piece_ids=tumor_piece_ids[!(tumor_piece_ids %in% c('',''))]
cell_types=paste(tumor_piece_ids,"Tumor",sep='_')
###Check that cell type distribution looks OK:
table(ATAC@meta.data %>% select ('Piece_ID','cell_type'))


#BAP1-mutants vs non-mutants

DefaultAssay(ATAC) <- 'peaksMACS2'
Idents(ATAC)=ATAC$test
both_m=c('','')

all_da_peaks=NULL
for (cell_type in bap1_s){
da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_type,
  ident.2='NOT_BAP1_mutant',
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
write.table(all_da_peaks, paste("out/da_peaks_BAP1mut_vs_NonMutants.min.pct0.1.",
"20210713.tsv",sep=""),sep="\t",quote=FALSE,
row.names=FALSE)


