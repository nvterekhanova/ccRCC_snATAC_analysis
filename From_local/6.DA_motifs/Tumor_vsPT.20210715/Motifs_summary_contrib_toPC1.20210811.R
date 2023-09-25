library(Signac)
library(ChIPseeker)
library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(reshape)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
ecdf_fun = function(x,perc) ecdf(x)(perc)

setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/6.DA_motifs/Tumor_vsPT.20210715/')


#2021-05-12: now make separate NAT-samples in the plot:
da_peaks=read.table('../../14.Cell_of_origin/20210805_Correlation_heatmap/out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
da_peaks$Sample=gsub('(.*)_.*','\\1',da_peaks$cell_type)
da_peaks$cell_typ1=gsub('(.*)_(.*)','\\2',da_peaks$cell_type)
da_peaks=da_peaks[(da_peaks$Sample=='PT' | da_peaks$cell_typ1 %in% c('Tumor','EMT tumor cells')) & !is.na(da_peaks$cell_type),]
da_peaks=da_peaks[!(da_peaks$cell_typ1 %in% c(12,14)),]
both=dcast(da_peaks, TF_Name~cell_type,value.var='mean_score')
#both=both[,1:29]
rownames(both)=both[,1]
both=both[,-1]

tfs_pc_1=c('HNF1A','HNF1B','NR1I3','SRF','XBP1','HNF4G','HNF4A','TCF7L1','HNF4A(var.2)','HMBOX1','ZNF410','TCF7','BARX2','TCF7L2',
'NR1H2::RXRA','RXRG','RXRB','PPARD','TFAP4(var.2)','LEF1')[1:10]
n=24
to_plot=both[rownames(both) %in% tfs_pc_1,]

annot_col=as.data.frame(colnames(both))
colnames(annot_col)=('ID')
annot_col$ID_1=gsub('(.*)_(.*)','\\1',annot_col$ID)
annot_col$ID_2=gsub('(.*)_(.*)','\\2',annot_col$ID)
annot_col$Cell_type=ifelse(annot_col$ID_1=='PT','PT',annot_col$ID_2)
annot_col$Sample=ifelse(annot_col$ID_1=='PT',paste("C",annot_col$ID_2,sep=''),annot_col$ID_1)

col_cell_t=c('PT'='#1B9E77','Tumor'='#E7298A','EMT tumor cells'='#E41A1C')
column_ha=HeatmapAnnotation(Cell_type=annot_col$Cell_type,IDs=anno_text(annot_col$Sample),col=list(Cell_type=col_cell_t))

x=Heatmap(to_plot,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score',show_column_dend=FALSE,show_row_dend=FALSE,bottom_annotation=column_ha, show_column_names=F)
#m=rowMeans(to_plot)
#to_plot_1=to_plot[names(m[m>0]),]
#x=Heatmap(to_plot_1,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score')


pdf("TopTFs_Scores_Separating_PT_Tumor_PC1.v2.20210902.pdf",width=10,height=4,useDingbats=F)
x
dev.off()


