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

setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/6.DA_motifs/For_cell_of_origin.20210721/Motifs_contributing_toPC_plot/')


#2021-05-12: now make separate NAT-samples in the plot:
da_peaks=read.table('../../../14.Cell_of_origin/20210805_Correlation_heatmap/out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
da_peaks$Sample=gsub('(.*)_.*','\\1',da_peaks$cell_type)
da_peaks$cell_typ1=gsub('(.*)_(.*)','\\2',da_peaks$cell_type)
da_peaks=da_peaks[(da_peaks$Sample=='PT' | da_peaks$cell_typ1 %in% c('Tumor','EMT tumor cells')) & !is.na(da_peaks$cell_type),]
da_peaks=da_peaks[!(da_peaks$cell_typ1 %in% c(12,14)),]
both=dcast(da_peaks, TF_Name~cell_type,value.var='mean_score')
#both=both[,1:29]
rownames(both)=both[,1]
both=both[,-1]

tfs_pc_1=c('HNF1A','HNF1B','NR1I3','SRF','XBP1','HNF4G','HNF4A','TCF7L1','HNF4A(var.2)','HMBOX1','ZNF410','TCF7','BARX2','TCF7L2',
'NR1H2::RXRA','RXRG','RXRB','PPARD','TFAP4(var.2)','LEF1')
n=24
to_plot=both[rownames(both) %in% tfs_pc_1,]
x=Heatmap(to_plot,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score')
#m=rowMeans(to_plot)
#to_plot_1=to_plot[names(m[m>0]),]
#x=Heatmap(to_plot_1,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score')


pdf("TopTFs_Scores_Separating_PT_Tumor_PC1.20210811.pdf",width=10,height=6,useDingbats=F)
x
dev.off()

######Now, try plotting ccRCC-tumor TFs
da_peaks=read.table('../../../14.Cell_of_origin/20210805_Correlation_heatmap/out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.20210804.tsv',sep='\t',header=TRUE)
da_peaks$Sample=gsub('(.*)_.*','\\1',da_peaks$cell_type)
da_peaks$cell_typ1=gsub('(.*)_(.*)','\\2',da_peaks$cell_type)
da_peaks=da_peaks[(da_peaks$Sample=='PT' | da_peaks$cell_typ1 %in% c('Tumor')) & !is.na(da_peaks$cell_type),]
da_peaks=da_peaks[!(da_peaks$cell_typ1 %in% c(12,14)),]
both=dcast(da_peaks, TF_Name~cell_type,value.var='mean_score')
#both=both[,1:29]
rownames(both)=both[,1]
both=both[,-1]


tfs_tumor=c('RELA','NFKB2','NFKB1','REL','HIF1A','KLF9','RREB1','ARNT::HIF1A','RBPJ','MXI1','ZNF75D','HSF2','SREBF2','NEUROD1','NEUROG2(var.2)','TBXT')

to_plot=both[rownames(both) %in% tfs_tumor,]
x=Heatmap(to_plot,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score')

pdf("TopTFs_Scores_Separating_PT_Tumor_PC1.20210811.pdf",width=10,height=6,useDingbats=F)
x
dev.off()




