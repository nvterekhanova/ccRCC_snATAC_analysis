####2021-08-05:
####2021-08-08, adding S1/2/3 annotation from Ruiyang-analysis:

library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(reshape2)
library(ggfortify)
library(tidyverse)


theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/14.Cell_of_origin/20210805_Correlation_heatmap/')



##########################################################################################
#######Now merging with the Average motif-activity score for each motif across cell types:
##########################################################################################

score=read.table('out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1=score_1[!(score_1$cell_type1 %in% c(12,14)),]
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:14),"PT",score_1$cell_type1)
score_1=score_1[!is.na(score_1$cell_type1),]

#score_1=score_1[score_1$cell_type1 %in% c('Tumor','Microglia','Oligodendrocytes','Fibroblasts','Neurons','T cells'),]

score=score_1
score=score[!(score$cell_type1 %in% c("Macrophages","NK cells","Unknown","CD8+ T-cells","CD4+ T-cells","DC","Endothelial cells","Fibroblasts")),]
score$cell_type1=ifelse(score$cell_type1=="Microglia","TAM",score$cell_type1)
score$cell_type=paste(score$Piece_ID,score$cell_type1,sep="_")
all_2=dcast(score,cell_type~TF_Name,value.var="mean_score")
rownames(all_2)=all_2[,1]
all_2=all_2[,-1]

mydata.cor = cor(t(all_2), method = c("spearman"))
#mydata.cor = cor(t(all_2), method = c("pearson"))

annot=score
annot=annot[!duplicated(annot$cell_type),]
rownames(annot)=annot$cell_type
annot=annot[rownames(mydata.cor),]

col_cell_t=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
names(col_cell_t)=c("Tumor","PT","Distal convoluted tubule","Intercalated cells","Loop of Henle","EMT tumor cells","Principle cells","Podocytes")
row_ha=rowAnnotation(Cell_type=annot$cell_type1,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t))


h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F)
pdf('Heatmap_Correlation.NAT_cellTypes.PTbyCluster.2021-08-05.pdf',width=12,height=10)
h
dev.off()

####Try only PT and Tumor cells:
score=read.table('out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1=score_1[!(score_1$cell_type1 %in% c(12,14)),]
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:14),"PT",score_1$cell_type1)
score_1=score_1[!is.na(score_1$cell_type1),]

#score_1=score_1[score_1$cell_type1 %in% c('Tumor','Microglia','Oligodendrocytes','Fibroblasts','Neurons','T cells'),]

score=score_1
score=score[score$cell_type1 %in% c('Tumor','PT'),]
score$cell_type=paste(score$Piece_ID,score$cell_type1,sep="_")
all_2=dcast(score,cell_type~TF_Name,value.var="mean_score")
rownames(all_2)=all_2[,1]
all_2=all_2[,-1]

mydata.cor = cor(t(all_2), method = c("spearman"))
#mydata.cor = cor(t(all_2), method = c("pearson"))

annot=score
annot=annot[!duplicated(annot$cell_type),]
rownames(annot)=annot$cell_type
annot=annot[rownames(mydata.cor),]

col_cell_t=c(brewer.pal(n = (length(unique(score$cell_type1))+1), name = "Paired"))
names(col_cell_t)=c(unique(score$cell_type1),'Other')
row_ha=rowAnnotation(Cell_type=annot$cell_type1,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t))

h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F)
pdf('Heatmap_Correlation.PT_Tumor.PTbyCluster.2021-08-05.pdf',width=12,height=10)
h
dev.off()


####Try making PCA-plot:
###############
#######HERE!###
###############
score=read.table('out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1=score_1[!(score_1$cell_type1 %in% c(12,14)),]
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:14),"PT",score_1$cell_type1)
score_1=score_1[!is.na(score_1$cell_type1),]

score=score_1
score=score[score$cell_type1 %in% c("PT","Tumor"),]
score$cell_type=paste(score$Piece_ID,score$cell_type1,sep="_")
all_2=dcast(score,cell_type~TF_Name,value.var="mean_score")
rownames(all_2)=all_2[,1]
all_2=all_2[,-1]

###make labels consistent with the main figg heatmap:
all_2$Cell_type=gsub('(.*)_(.*)','\\2',rownames(all_2))
all_2$Sample=gsub('(.*)_(.*)','\\1',rownames(all_2))
all_2$Sample=ifelse(all_2$Cell_type=='PT',paste('C',all_2$Sample,sep=''),'')
pca_res <- prcomp(all_2[1:(ncol(all_2)-2)], scale. = TRUE)
autoplot(pca_res,data=all_2,colour='Cell_type',label=T,label.label='Sample')

pdf('PCA_plot_ccRCC_Tumor_PTbyClusster.v2.20210902.pdf',width=5.5,height=4.5)
autoplot(pca_res,data=all_2,colour='Cell_type',label=T,label.label='Sample',label.vjust=1.5,shape='Cell_type')+theme_classic()+scale_color_manual(values=c('blue','red'))+scale_shape_manual(values=c(19,15))
#+scale_color_manual(values=col_cell_t)
dev.off()

####Identify TFs contributing to Dims 1/2:
library("factoextra")
#axes corresponds to dim:
dim_1=fviz_contrib(pca_res, choice="var", axes = 1, top=20)
pdf('Contribution_toPC1_ccRCC_Tumor_PTbyCluster_consistent_with_snRNA.v2.20210902.pdf',width=5.8,height=2.6)
print(dim_1)
dev.off()
dim_2=fviz_contrib(pca_res, choice="var", axes = 2, top=20)
pdf('Contribution_toPC2_ccRCC_Tumor_PTbyCluster_consistent_with_snRNA.v2.20210902.pdf',width=6.12,height=2.6)
print(dim_2)
dev.off()



#####Also try adding fibroblasts:
####Try only PT and Tumor cells:
score=read.table('out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1=score_1[!(score_1$cell_type1 %in% c(12,14)),]
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:14),"PT",score_1$cell_type1)
score_1=score_1[!is.na(score_1$cell_type1),]


score=score_1
score=score[score$cell_type1 %in% c("PT","Tumor","Fibroblasts","EMT tumor cells"),]
score$cell_type1=ifelse(score$cell_type1=="Microglia","TAM",score$cell_type1)
score$cell_type=paste(score$Piece_ID,score$cell_type1,sep="_")
all_2=dcast(score,cell_type~TF_Name,value.var="mean_score")
rownames(all_2)=all_2[,1]
all_2=all_2[,-1]

mydata.cor = cor(t(all_2), method = c("spearman"))
#mydata.cor = cor(t(all_2), method = c("pearson"))

score$Sample=gsub('(.*-.*)-.*','\\1',score$Piece_ID)
clin=read_delim('../../clinical_Pan-cancer.July2021.tsv',delim='\t')
clin=as.data.frame(clin)
clin=clin[clin[,1]=='CCRCC' & clin$case_id %in% score$Sample,]
clin=clin[,c(2,32,44)]
clin$Grade=gsub('(.*):.*','\\1',clin[,2])
colnames(clin)[1]='Sample'
annot=score
annot=annot[!duplicated(annot$cell_type),]
annot=merge(annot,clin,all.x=T)
rownames(annot)=annot$cell_type
annot=annot[rownames(mydata.cor),]
colnames(annot)[8]='Stage'
annot$Stage=ifelse(annot$cell_type1 %in% c('Fibroblasts','EMT Tumor cells'),NA,annot$Stage)
annot$Grade=ifelse(annot$cell_type1 %in% c('Fibroblasts','EMT Tumor cells'),NA,annot$Grade)
annot$Segment=NA
annot$Segment=ifelse(annot$cell_type1=='PT' & annot$Sample %in% c(0,5,8),"S1/S2",annot$Segment)
annot$Segment=ifelse(annot$cell_type1=='PT' & annot$Sample %in% c(13),"S2/S3",annot$Segment)
annot$Segment=ifelse(annot$cell_type1=='PT' & annot$Sample %in% c(3,4,11),"No_segm_markers",annot$Segment)
segm_cols=c('#005824','#41ae76','#edf8fb')
names(segm_cols)=c('S1/S2','S2/S3','No_segm_markers')

col_cell_t=c(brewer.pal(n = (length(unique(score$cell_type1))-1), name = "Paired"),"red")
names(col_cell_t)=c(unique(score$cell_type1)[c(1,3,2,4)])
cols_grade=c(brewer.pal(n=4,name='YlGnBu'))
names(cols_grade)=c('G1','G2','G3','G4')
cols_stage=c(brewer.pal(n=4,name='YlOrRd'))
names(cols_stage)=c('Stage I','Stage II','Stage III','Stage IV')
row_ha=rowAnnotation(Cell_type=annot$cell_type1,PT_segment=annot$Segment,Hist_grade=annot$Grade,Stage=annot$Stage,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t,Hist_grade=cols_grade,
Stage=cols_stage, PT_segment=segm_cols))

h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F)
pdf('Heatmap_Correlation.PT_Tumor.EMT_Fibroblasts_PTbyCluster.Grade.Stageadded.PTSegmAnnotation.2021-08-05.pdf',width=12,height=10)
h
dev.off()


