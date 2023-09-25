#####2022-02-13
library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(reshape2)
library(ChIPseeker)
library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Signac)


setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/5.DA_peaks/2021-07-12_BAP1_PBRM1_vs_Nonmutants_28samples/BAP1_peaks_for_ImmuneStroma.revision/')
#tab_1=read.table('../BAP1_vs_NonMutants/out/DOWN_BAP1mutants_vs_nonMutans.Accessibility.20210713.tsv',sep='\t',header=TRUE)
tab_1=read.table('../BAP1_vs_NonMutants/counts.Norm.forPlotting.20211011/DOWN_BAP1mutants_vs_nonMutans.Accessibility.20211011.tsv',sep='\t',header=TRUE)

row.names(tab_1)=paste(tab_1[,1],tab_1[,2],sep='_')
#tab_2=read.table('../BAP1_vs_NonMutants/out/UP_BAP1mutants_vs_nonMutans.Accessibility.20210713.tsv',sep='\t',header=TRUE)
tab_2=read.table('../BAP1_vs_NonMutants/counts.Norm.forPlotting.20211011/UP_BAP1mutants_vs_nonMutans.Accessibility.20211011.tsv',sep='\t',header=TRUE)
row.names(tab_2)=paste(tab_2[,1],tab_2[,2],sep='_')
tab_2=tab_2[rownames(tab_1),]
colnames(tab_1)=gsub('\\.','\\-',colnames(tab_1))
colnames(tab_2)=gsub('\\.','\\-',colnames(tab_2))

tab_2=tab_2[,3:ncol(tab_2)]
tab=cbind(tab_1,tab_2)

tab$Group.1=gsub('_','-',tab$Group.1)
tab$Group.2=gsub('_','.',tab$Group.2)
colnames(tab)=gsub('\\.','\\-',colnames(tab))
tab$ID=paste(tab[,1],tab[,2],sep='_')
rownames(tab)=tab$ID
colnames(tab)[1:2]=c('Case','Cell_type')
#tab=tab[tab$Cell_type %in% c('Tumor','PT'),]

#tab=tab[tab$Cell_type %in% c('Cell.line','Tumor','PT','Immune'),]

tab=tab[,3:(ncol(tab)-1)]

annot=read.table('../../2021-06-10_BAP1_specific/Accessibility/Sample_categories.20210503.txt',sep='\t',header=TRUE)

annot=as.data.frame(rownames(tab))
colnames(annot)='ID'
rownames(annot)=annot$ID
annot$Case=gsub('(.*)_(.*)','\\1',rownames(annot))
annot$Cell_type=gsub('(.*)_(.*)','\\2',rownames(annot))
annot$BAP1_mutation=ifelse(annot$Case %in% c(bap1_s,'BAP1-786O'),'BAP1_mutant','NOT_BAP1_mutant')
annot$PBRM1_mutation=ifelse(annot$Case %in% pbrm1_s,'PBRM1_mutant','NOT_PBRM1_mutant')

#rownames(tab)=gsub('(.*)_.*','\\1',rownames(tab))



###Try remove cell types with low cell count:
meta=read.table('28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=T)
meta$cell_group=meta$cell_type
meta$cell_group=ifelse(meta$cell_group %in% c('DC','CD4+ T-cells','Macrophages','NK cells','CD8+ T-cells'),'Immune',meta$cell_group)
meta$cell_group=ifelse(meta$cell_group %in% c('Endothelial cells','Fibroblasts'),'Stroma',meta$cell_group)
meta_s=meta %>% dplyr::select ('cell_group','Piece_ID')
meta_s$Count=1
meta_s1=aggregate(meta_s$Count,by=list(meta_s$Piece_ID,meta_s$cell_group),FUN='sum')
colnames(meta_s1)=c('Piece_ID','Cell_group','Cell_count')
meta_s1=meta_s1[meta_s1$Cell_group %in% annot$Cell_type,]
meta_s1$ID=paste(meta_s1$Piece_ID,meta_s1$Cell_group,sep='_')
meta_s1=meta_s1[meta_s1$Cell_count>100,]

annot=annot[annot$ID %in% meta_s1$ID,]
tab=tab[rownames(annot),]
row_ha=rowAnnotation(Cell_type=annot$Cell_type, BAP1_status=annot$BAP1_mutation,PBRM1_status=annot$PBRM1_mutation,ID=anno_text(annot$Case),col=list(BAP1_status=c('BAP1_mutant'='red','NOT_BAP1_mutant'='white smoke'),PBRM1_status=c('PBRM1_mutant'='blue','NOT_PBRM1_mutant'='white smoke'),Cell_type=c('Cell.line'='#a65628','Immune'='#ff7f00','PT'='#4daf4a','Tumor'='#984ea3','Stroma'='#ffff33')))

col_annot=as.data.frame(colnames(tab))
colnames(col_annot)='peak'
col_annot$Type=ifelse(col_annot$peak %in% colnames(tab_1),'DOWN','UP')
x=Heatmap(scale(tab),show_row_names = F,show_column_names = FALSE,show_row_dend=FALSE,show_column_dend=FALSE,right_annotation=row_ha,name='Peak_accessibility', row_split = factor(annot$Cell_type),column_split=factor(col_annot$Type))

#cluster_row_slices=F, column_split = factor(c(rep("DOWN", 3890),rep("UP",647)),levels=c('UP','DOWN')),cluster_column_slices=F,use_raster=FALSE)

pdf("UP_and_DOWN_BAP1_specific_peaks.countsNorm.Incl.ImmuneStroma.More100cells.20220213.pdf",width=13,height=10)
print(x)
dev.off()
