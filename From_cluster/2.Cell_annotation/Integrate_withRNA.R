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
library(plyr)
library(RColorBrewer)

library(EnsDb.Hsapiens.v86)
args = commandArgs(trailingOnly=TRUE)

sample=args[1]

print(sample)


date=args[2]
dir.create(paste("out/",sample,sep=''))
dir.create(paste("out/",sample,'/Annot_',date,sep=''))

RNA=readRDS(paste('../../snRNA_processed_Yige/snRNA_inputs/',sample,'_processed.rds',sep=''))
ATAC=readRDS(paste("../1.Create_rds/out/",sample,"/",sample,"_processed_atac.rds",sep=""))
DefaultAssay(ATAC)='X500peaksMACS2'

#####Labelling cell-types for the snRNA

meta=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/',
'snRNA_processed_Yige/meta_data.20200505.v1.tsv',sep=''),sep='\t',header=TRUE)
case=unique(as.character((meta$Aliquot.snRNA.WU[meta$Aliquot.snRNA==sample]))) 

#annot=read.table('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/snRNA_processed_Yige/31Aliquot.Barcode2CellType.20201121.v1.tsv',header=TRUE,sep="\t")

annot=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/',
'snRNA_processed_Yige/33Aliquot.Barcode2CellType.20210423.v1.tsv',sep=''),header=TRUE,sep="\t")


cols=c(brewer.pal(n = 12, name = "Paired"),"grey","#8B008B","#660000")
names(cols)= unique(as.character(annot$Cell_group14_w_transitional))

annot=annot[annot$orig.ident==sample,]

orig_1=RNA$orig.ident
orig_1=as.data.frame(orig_1)

#orig_1$cell_id=rownames(orig_1)
orig_1$individual_barcode=rownames(orig_1)


res=merge(orig_1, annot,all.x=TRUE)

#RNA$cell_type=res$cell_type
#RNA$cell_type=res$Cell_type.shorter
RNA$cell_type=res$Cell_group14_w_transitional



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
  refdata = RNA$cell_type,
  weight.reduction = ATAC[['lsi']],
  dims=2:30
)

ATAC <- AddMetaData(object = ATAC, metadata = predicted.labels)
ATAC$Piece_ID=case
###Making plots:

plot1 <- DimPlot(
  object = RNA,
  group.by = 'cell_type',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(case,' snRNA-seq',sep=""))+scale_color_manual(values=cols)

plot2 <- DimPlot(
  object = ATAC,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(case,' snATAC-seq',sep=""))+scale_color_manual(values=cols)

pdf(paste("out/",sample,"/Annot_",date,"/",sample,"_",case,"_snRNA_snATAC_integrated.pdf",sep=""),height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()

atac_meta=ATAC@meta.data
write.table(atac_meta,paste("out/",sample,"/Annot_",date,"/",sample,"_cellTyped.meta.data",sep=""),sep='\t',
quote=FALSE)
#saveRDS(ATAC,paste("out/",sample, "/Annot_",date,"/", sample,"_ATAC_cellTyped.rds",sep=""))
