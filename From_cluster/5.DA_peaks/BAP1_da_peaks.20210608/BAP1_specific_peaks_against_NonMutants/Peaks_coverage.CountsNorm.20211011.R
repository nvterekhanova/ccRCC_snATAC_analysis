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



ATAC=readRDS(paste('../../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'20210707.rds',sep=''))


###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../../../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type


####First DOWN-peaks:
peaks=read_delim('data/DOWN_BAP1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210713.tsv',
delim='\t')
peaks=as.data.frame(peaks)

#ATAC=atac
x=ATAC@assays$peaksMACS2@counts
x1=ATAC@assays$peaksMACS2@counts
x1=x1[rownames(x1) %in% peaks$peak,]
x1=as.data.frame(x1)
y=t(x1)
y=as.data.frame(y)

#ATAC$cell_type=ifelse(ATAC$Piece_ID %in% c('BAP1_786O','Control_786O'),'cell_line',ATAC$cell_type)
y$Case=ATAC$Piece_ID
y$Cell_type=ATAC$cell_type

m=colSums(x)
y1=y
y1[,1:(ncol(y1)-2)]=y1[,1:(ncol(y1)-2)]/m
y1[,1:(ncol(y1)-2)]=y1[,1:(ncol(y1)-2)]*10000

y=y1
y$Cell_group='Other'
y$Cell_group=ifelse(y$Cell_type %in% c('DC','CD4+ T-cells','CD8+ T-cells','Macrophages','NK cells'),
'Immune',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type	%in% c('Endothelial cells','Fibroblasts'),'Stroma',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='Tumor','Tumor',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='PT','PT',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='Loop of Henle','LOH',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='cell_line','Cell_line',y$Cell_group)
y=y[y$Cell_group!='Other',]


###Now aggregate:
t1=aggregate(y[,1:(ncol(y)-3)],by=list(y$Case, y$Cell_group),FUN='mean',na.action = na.omit)

write.table(t1,"out/Accessibility.Counts/DOWN_BAP1mutants_vs_nonMutans.Accessibility.20211011.tsv",sep='\t',
quote=F,row.names=F)

----------------------------------
#Now UP peaks:
----------------------------------
peaks=read_delim('data/UP_BAP1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210713.tsv',
delim='\t')
peaks=as.data.frame(peaks)

#ATAC=atac
x=ATAC@assays$peaksMACS2@counts

x1=ATAC@assays$peaksMACS2@counts
x1=x1[rownames(x1) %in% peaks$peak,]
x1=as.data.frame(x1)
y=t(x1)
y=as.data.frame(y)

#ATAC$cell_type=ifelse(ATAC$Piece_ID %in% c('BAP1_786O','Control_786O'),'cell_line',ATAC$cell_type)
y$Case=ATAC$Piece_ID
y$Cell_type=ATAC$cell_type

m=colSums(x)
y1=y
y1[,1:(ncol(y1)-2)]=y1[,1:(ncol(y1)-2)]/m
y1[,1:(ncol(y1)-2)]=y1[,1:(ncol(y1)-2)]*10000

y=y1
y$Cell_group='Other'
y$Cell_group=ifelse(y$Cell_type %in% c('DC','CD4+ T-cells','CD8+ T-cells','Macrophages','NK cells'),
'Immune',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type	%in% c('Endothelial cells','Fibroblasts'),'Stroma',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='Tumor','Tumor',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='PT','PT',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='Loop of Henle','LOH',y$Cell_group)
y$Cell_group=ifelse(y$Cell_type=='cell_line','Cell_line',y$Cell_group)
y=y[y$Cell_group!='Other',]


###Now aggregate:
t1=aggregate(y[,1:(ncol(y)-3)],by=list(y$Case, y$Cell_group),FUN='mean',na.action = na.omit)

write.table(t1,"out/Accessibility.Counts/UP_BAP1mutants_vs_nonMutans.Accessibility.20211011.tsv",sep='\t',quote=F)
