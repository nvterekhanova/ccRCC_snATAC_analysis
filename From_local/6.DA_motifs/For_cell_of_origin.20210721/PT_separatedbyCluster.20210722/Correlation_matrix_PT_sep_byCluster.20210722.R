####2021-04-18:
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

setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/6.DA_motifs/For_cell_of_origin.20210721/PT_separatedbyCluster.20210722/')



##########################################################################################
#######Now merging with the Average motif-activity score for each motif across cell types:
##########################################################################################

###############
#######HERE!###
###############
score=read.table('Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210721.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:11),"PT",score_1$cell_type1)

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

col_cell_t=c(brewer.pal(n = length(unique(score$cell_type1)), name = "Paired"))
names(col_cell_t)=unique(score$cell_type1)
row_ha=rowAnnotation(Cell_type=annot$cell_type1,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t))


h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F)
pdf('Heatmap_Correlation.NAT_cellTypes.PTbyCluster.2021-07-22.pdf',width=12,height=10)
h
dev.off()

####Try only PT and Tumor cells:
score=read.table('Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210721.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:11),"PT",score_1$cell_type1)


score=score_1
score=score[score$cell_type1 %in% c("PT","Tumor","Fibroblast"),]
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

col_cell_t=c(brewer.pal(n = length(unique(score$cell_type1)), name = "Paired"))
names(col_cell_t)=c(unique(score$cell_type1),"Other")
row_ha=rowAnnotation(Cell_type=annot$cell_type1,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t))

h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F)
pdf('Heatmap_Correlation.PT_Tumor.PTbyCluster.2021-07-22.pdf',width=12,height=10)
h
dev.off()


###########################
####Try making PCA-plot:###
###########################

score=read.table('Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210721.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:11),"PT",score_1$cell_type1)


score=score_1
score=score[score$cell_type1 %in% c("PT","Tumor"),]
score$cell_type1=ifelse(score$cell_type1=="Microglia","TAM",score$cell_type1)
score$cell_type=paste(score$Piece_ID,score$cell_type1,sep="_")
all_2=dcast(score,cell_type~TF_Name,value.var="mean_score")
rownames(all_2)=all_2[,1]
all_2=all_2[,-1]

all_2$Cell_type=gsub('(.*)_(.*)','\\2',rownames(all_2))
all_2$Sample=gsub('(.*)_(.*)','\\1',rownames(all_2))
pca_res <- prcomp(all_2[1:(ncol(all_2)-2)], scale. = TRUE)
autoplot(pca_res,data=all_2,colour='Cell_type',label=T,label.label='Sample')

pdf('PCA_plot_ccRCC_Tumor_PTbyClusster.20210722.pdf',width=9,height=8)
autoplot(pca_res,data=all_2,colour='Cell_type',label=T,label.label='Sample',label.vjust=1.5)+theme_classic()+scale_color_manual(values=c('blue','red'))
#+scale_color_manual(values=col_cell_t)
dev.off()

#####Also try adding fibroblasts:
####Try only PT and Tumor cells:
setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/6.DA_motifs/For_cell_of_origin.20210721/PT_separatedbyCluster.20210722/')
score=read.table('Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210721.tsv',sep='\t',header=TRUE)
score_1=score
score_1$Piece_ID=gsub('(.*)_.*','\\1',score_1$cell_type)
score_1$cell_type1=gsub('(.*)_(.*)','\\2',score_1$cell_type)
score_1$Piece_ID=ifelse(score_1$Piece_ID=='PT',score_1$cell_type1,score_1$Piece_ID)
score_1$cell_type1=ifelse(score_1$cell_type1 %in% c(0:11),"PT",score_1$cell_type1)


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
clin=read_delim('../../../clinical_Pan-cancer.July2021.tsv',delim='\t')
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

col_cell_t=c(brewer.pal(n = (length(unique(score$cell_type1))-1), name = "Paired"),"red")
names(col_cell_t)=c(unique(score$cell_type1)[c(1,3,2,4)])
cols_grade=c(brewer.pal(n=4,name='YlGnBu'))
names(cols_grade)=c('G1','G2','G3','G4')
cols_stage=c(brewer.pal(n=4,name='YlOrRd'))
names(cols_stage)=c('Stage I','Stage II','Stage III','Stage IV')
row_ha=rowAnnotation(Cell_type=annot$cell_type1,Hist_grade=annot$Grade,Stage=annot$Stage,IDs=anno_text(annot$Piece_ID),col=list(Cell_type=col_cell_t,Hist_grade=cols_grade,
Stage=cols_stage))

h=Heatmap(mydata.cor,col= colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),show_column_names = FALSE,name="Motif score",right_annotation=row_ha,show_row_names=F)
pdf('Heatmap_Correlation.PT_Tumor.EMT_Fibroblasts_PTbyCluster.Grade.Stageadded.2021-07-22.pdf',width=12,height=10)
h
dev.off()


#####Now try to make plots of Grade across different groups
tab=annot[annot$cell_type1=='Tumor',]
mesenchymal=c('')
close_pt=c('')
close_pt_2=c('')

tab$Tumor_group='Intermediate'
tab$Tumor_group=ifelse(tab$Piece_ID %in% close_pt,'Close_toPT',tab$Tumor_group)
tab$Tumor_group=ifelse(tab$Piece_ID %in% close_pt_2,'Close_toPT_group2',tab$Tumor_group)
tab$Tumor_group=ifelse(tab$Piece_ID %in% mesenchymal,'Mesenchymal-like',tab$Tumor_group)


tab$Hist_Grade=as.numeric(gsub('G(.*)','\\1',tab$Grade))

p <- ggplot(data = tab, aes(x=Tumor_group,y=Hist_Grade)) +

p <- p + geom_boxplot(outlier.shape = NA,width=0.6) 

p <- p + geom_jitter(shape=19, position=position_jitter(0.2),size=1,aes(color=Grade),alpha=0.5)
             
p <- p + theme_bw() + theme_nogrid() + labs(title="",x="Tumor sample group",y="Histologic grade")

p <- p + theme(axis.text.x = element_text(colour="black", size=10, ), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())

p <- p + theme(legend.position = "none") + scale_color_manual(values=cols)
