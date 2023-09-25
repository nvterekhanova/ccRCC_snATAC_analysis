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
#annot=read.table('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/Merged_snATAC/Merge.20210503_manualReannotation/26_ccRCC_snATAC_ManualReviwed.v2.20210509.tsv',sep='\t',header=TRUE)
da_peaks=read.table('out/Score_difference.Tumor_Normal_comparison.20210715.tsv',sep='\t',header=TRUE)
up=da_peaks[da_peaks$diff>0 & da_peaks$FDR<0.05,]
up=up %>% dplyr::select ('cell_t1','TF_Name','diff')
colnames(up)[1]='Sample'
up=dcast(up,TF_Name~Sample,value.var='diff')
rownames(up)=up[,1]
up=up[,-1]
up_1 =up %>% mutate(Count_up = rowSums(!is.na(up)))
colnames(up_1)[1:24]=paste(colnames(up_1)[1:24],'Diff_vs_PT',sep='_')
up_1$Mean_diff=rowMeans(up_1[,1:24],na.rm=TRUE)
up_1=up_1[order(-up_1$Mean_diff),]
up_1=up_1[order(-up_1$Count_up),]

write.table(up_1,'out/Motif_scores_Tumor_samples_vs_PT.20210715.tsv',sep='\t',quote=FALSE)


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
to_plot=both[rownames(both) %in% rownames(up_1[up_1$Count_up>=n,]),]
x=Heatmap(to_plot,col= colorRamp2(c(-1.5, 0, 1.5), c("#377EB8", "white", "#E41A1C")),row_order=rownames(up_1[up_1$Count_up>=n,]),name='Motif score')
#x=Heatmap(to_plot,col= colorRamp2(c(-1.5, 0, 1.5), c("white", "#ffeda0","#bd0026")),row_order=rownames(up_1[up_1$Count_up>=n,]))

#pdf("TopTFs_Scores_Tumor_PT.20210510.pdf",width=10,height=7,useDingbats=F)
pdf("TopTFs_Scores_Tumor_PT.N_24.20210715.pdf",width=10,height=6,useDingbats=F)
x
dev.off()

###Plot difference only:
all=da_peaks
tumor=all %>% select ('cell_t1','TF_Name','mean_score1','diff')
colnames(tumor)[c(1,3)]=c('Sample','mean_score')
both=tumor

both=dcast(both,TF_Name~Sample,value.var='diff')
rownames(both)=both[,1]
both=both[,-1]
#both$PT=0

n=23
to_plot=both[rownames(both) %in% rownames(up_1[up_1$Count_up>=n,]),]
x=Heatmap(to_plot,col= colorRamp2(c(-1, 0, 2), c("#377EB8", "white", "#E41A1C")),row_order=rownames(up_1[up_1$Count_up>=n,]),name='Score_Difference_vs_PT')
pdf("TopTFs_score_difference_fromPT.20210510.pdf",width=10,height=7,useDingbats=F)
x
dev.off()


####Try n=24:
n=24
to_plot=both[rownames(both) %in% rownames(up_1[up_1$Count_up>=n,]),]
x=Heatmap(to_plot,col= colorRamp2(c(-1, 0, 2), c("#377EB8", "white", "#E41A1C")),row_order=rownames(up_1[up_1$Count_up>=n,]),name='Score_Difference_vs_PT')
pdf("TopTFs_score_difference_fromPT.N_24.20210510.pdf",width=10,height=6,useDingbats=F)
x
dev.off()

########################
####ENDS HERE FOR NOW###
########################

#2021-05-12: now make separate NAT-samples in the plot:
da_peaks=read.table('out/Motif_score_perCell_group.20210715.tsv',sep='\t',header=TRUE)
both=dcast(da_peaks, TF_Name~cell_type,value.var='mean_score')
both=both[,1:29]
rownames(both)=both[,1]
both=both[,-1]

n=24
to_plot=both[rownames(both) %in% rownames(up_1[up_1$Count_up>=n,]),]
x=Heatmap(to_plot,col= colorRamp2(c(-1.5, 0, 1.5), c("#377EB8", "white", "#E41A1C")),row_order=rownames(up_1[up_1$Count_up>=n,]),name='Motif score')
#x=Heatmap(to_plot,col= colorRamp2(c(-1.5, 0, 1.5), c("white", "#ffeda0","#bd0026")),row_order=rownames(up_1[up_1$Count_up>=n,]))

#pdf("TopTFs_Scores_Tumor_PT.20210510.pdf",width=10,height=7,useDingbats=F)
pdf("TopTFs_Scores_Tumor_PT.N_24_2NAT.20210715.pdf",width=10,height=6,useDingbats=F)
x
dev.off()

###Plot difference only:
all=da_peaks
tumor=all %>% select ('cell_t1','TF_Name','mean_score1','diff')
colnames(tumor)[c(1,3)]=c('Sample','mean_score')
both=tumor

both=dcast(both,TF_Name~Sample,value.var='diff')
rownames(both)=both[,1]
both=both[,-1]
#both$PT=0

n=23
to_plot=both[rownames(both) %in% rownames(up_1[up_1$Count_up>=n,]),]
x=Heatmap(to_plot,col= colorRamp2(c(-1, 0, 2), c("#377EB8", "white", "#E41A1C")),row_order=rownames(up_1[up_1$Count_up>=n,]),name='Score_Difference_vs_PT')
pdf("TopTFs_score_difference_fromPT.20210510.pdf",width=10,height=7,useDingbats=F)
x
dev.off()


####Try n=24:
n=24
to_plot=both[rownames(both) %in% rownames(up_1[up_1$Count_up>=n,]),]
x=Heatmap(to_plot,col= colorRamp2(c(-1, 0, 2), c("#377EB8", "white", "#E41A1C")),row_order=rownames(up_1[up_1$Count_up>=n,]),name='Score_Difference_vs_PT')
pdf("TopTFs_score_difference_fromPT.N_24.20210510.pdf",width=10,height=6,useDingbats=F)
x
dev.off()

##############################
###ENDS HERE FOR NOOW#########
##############################
