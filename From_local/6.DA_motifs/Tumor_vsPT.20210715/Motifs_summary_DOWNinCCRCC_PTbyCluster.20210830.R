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


####pt-specific tfs:
setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/6.DA_motifs/Tumor_vsPT.20210715/')
da_peaks=read.table('out/Score_difference.Tumor_Normal_comparison.20210715.tsv',sep='\t',header=TRUE)
up=da_peaks[da_peaks$diff<0 & da_peaks$FDR<0.05,]
up=up %>% dplyr::select ('cell_t1','TF_Name','diff')
colnames(up)[1]='Sample'
up=dcast(up,TF_Name~Sample,value.var='diff')
rownames(up)=up[,1]
up=up[,-1]
up_1 =up %>% mutate(Count_up = rowSums(!is.na(up)))
colnames(up_1)[1:24]=paste(colnames(up_1)[1:24],'Diff_vs_PT',sep='_')
up_1$Mean_diff=rowMeans(up_1[,1:24],na.rm=TRUE)
up_1=up_1[order(up_1$Mean_diff),]
up_1=up_1[order(-up_1$Count_up),]

write.table(up_1,'out/DOWN_Motif_scores_Tumor_samples_vs_PT.20210720.tsv',sep='\t',quote=FALSE)


####now make heatmap:
all=da_peaks
pt=all %>% dplyr::select ('TF_Name','mean_score2')
pt$Sample='PT'
pt=pt[!duplicated(pt),]
colnames(pt)[2]='mean_score'
tumor=all %>% dplyr::select ('cell_t1','TF_Name','mean_score1')
#tumor=all %>% select ('cell_t1','TF_Name','mean_score1','diff')
colnames(tumor)[c(1,3)]=c('Sample','mean_score')

both=rbind(tumor,pt)
both=dcast(both,TF_Name~Sample,value.var='mean_score')
rownames(both)=both[,1]
both=both[,-1]

n=24
pt_specific_tfs=rownames(up_1[up_1$Count_up>=n,])




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
'NR1H2::RXRA','RXRG','RXRB','PPARD','TFAP4(var.2)','LEF1')
n=24
to_plot=both[rownames(both) %in% pt_specific_tfs,]
x=Heatmap(to_plot,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score')
#m=rowMeans(to_plot)
#to_plot_1=to_plot[names(m[m>0]),]
#x=Heatmap(to_plot_1,col= colorRamp2(c(-4, 0, 4), c("#377EB8", "white", "#E41A1C")),name='Motif score')


pdf("TopTFs_PTspecific_PTbyCluster.20210830.pdf",width=10,height=6,useDingbats=F)
x
dev.off()


