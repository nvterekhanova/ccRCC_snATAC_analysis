library(Signac)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(reshape)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(karyoploteR)
library(tidyverse)
ecdf_fun = function(x,perc) ecdf(x)(perc)


setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/5.DA_peaks/2021-07-12_ccRCC_specific_28samples/')
date='20210811'

da_peaks=read_delim('da_peaks_vs_PT_NAT.minPct0.1.20210712.tsv',delim='\t')
da_peaks=as.data.frame(da_peaks)

da_peaks=da_peaks[da_peaks$p_val_adj<0.05,]
da_peaks=da_peaks %>% dplyr::select ('peak','Sample','avg_log2FC')
up=da_peaks[da_peaks$avg_log2FC>0,]
up=dcast(up,peak~Sample,value.var='avg_log2FC')
rownames(up)=up[,1]
up=up[,-1]
up_1 =up %>% mutate(Count_up = rowSums(!is.na(up)))
up_1=up_1[order(-up_1$Count_up),]
colnames(up_1)[1:(ncol(up_1)-1)]=paste(colnames(up_1)[1:(ncol(up_1)-1)],'Signif_avg_lnFC',sep='_')
up_1$peak=rownames(up_1)


da_peaks=read_delim('da_peaks_vs_PT_NAT.minPct0.1.20210712.tsv',delim='\t')
da_peaks=as.data.frame(da_peaks)
da_peaks=da_peaks %>% dplyr::select ('peak','Sample','avg_log2FC')
up=dcast(da_peaks,peak~Sample,value.var='avg_log2FC')
colnames(up)[2:(ncol(up))]=paste('NonSignif_avg_lnFC',colnames(up)[2:(ncol(up))],sep='_')

res=merge(up_1,up,all.x=TRUE)


###cnv_corrected
da_p_c=read_delim('DA_peaks_Tumor_vs_PT_correctedbyCNV.20210808.tsv',delim='\t')
da_p_c=as.data.frame(da_p_c)
colnames(da_p_c)[4]=paste('CNV_corr_',colnames(da_p_c)[4],sep='')
da_p_c=da_p_c %>% dplyr::select ('peak','CNV_corr_p_adjust_bonf')
res=merge(res,da_p_c,all.x=TRUE)
res=res[order(-res$Count_up),]

res_1=res[res$Count_up>=12,]
#nrow(res_1)=4,514

res_2=res_1
for (i in (24+3):(24+2+24)){
	res_2=res_2[res_2[,i]>0 | is.na(res_2[,i]),]
}
#nrow(res_2)=1,538

res_3=res_2[res_2$CNV_corr_p_adjust_bonf<0.05 | is.na(res_2$CNV_corr_p_adjust_bonf),]
nrow(res_3)
#nrow(res_3)=1,526
res_3$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_3$peak)

to_plot=table(res_3$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('UP_Tumor_specific_DA_peaks.After_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='UP peaks after CNV correction',ylim=c(0,250))
dev.off()

res_2$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_2$peak)
to_plot=table(res_2$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('UP_Tumor_specific_DA_peaks.BeforeCNV_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='UP peaks before CNV correction',ylim=c(0,250))
dev.off()

peaks_1=StringToGRanges(res_3$peak, sep = c("-", "-"))
 ###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks=res_3
peaks$Gene=anno$SYMBOL
peaks$Type=anno$annotation
peaks$geneId=anno$geneId
peaks$peak_distanceToTSS=anno$distanceToTSS

#nrow(peaks[peaks$Type=='Promoter',])=510 510/1.526=21.8%

table(peaks$chr_peak[peaks$Type=='Promoter'])

write.table(peaks, "out/UP_Tumor_vsPT.Filtered.CNV_corrected.Annotated.20210811.tsv",sep='\t',row.names=F,quote=F)

---------------------------------
###############
#NOW DOWN PEAKS
###############
---------------------------------
da_peaks=read_delim('da_peaks_vs_PT_NAT.minPct0.1.20210712.tsv',delim='\t')
da_peaks=as.data.frame(da_peaks)

da_peaks=da_peaks[da_peaks$p_val_adj<0.05,]
da_peaks=da_peaks %>% dplyr::select ('peak','Sample','avg_log2FC')
up=da_peaks[da_peaks$avg_log2FC<0,]
up=dcast(up,peak~Sample,value.var='avg_log2FC')
rownames(up)=up[,1]
up=up[,-1]
up_1 =up %>% mutate(Count_down = rowSums(!is.na(up)))
up_1=up_1[order(-up_1$Count_down),]
colnames(up_1)[1:(ncol(up_1)-1)]=paste(colnames(up_1)[1:(ncol(up_1)-1)],'Signif_avg_lnFC',sep='_')
up_1$peak=rownames(up_1)


da_peaks=read_delim('da_peaks_vs_PT_NAT.minPct0.1.20210712.tsv',delim='\t')
da_peaks=as.data.frame(da_peaks)
da_peaks=da_peaks %>% dplyr::select ('peak','Sample','avg_log2FC')
up=dcast(da_peaks,peak~Sample,value.var='avg_log2FC')
colnames(up)[2:(ncol(up))]=paste('NonSignif_avg_lnFC',colnames(up)[2:(ncol(up))],sep='_')

res=merge(up_1,up,all.x=TRUE)


###cnv_corrected
da_p_c=read_delim('DA_peaks_Tumor_vs_PT_correctedbyCNV.20210808.tsv',delim='\t')
da_p_c=as.data.frame(da_p_c)
colnames(da_p_c)[4]=paste('CNV_corr_',colnames(da_p_c)[4],sep='')
da_p_c=da_p_c %>% dplyr::select ('peak','CNV_corr_p_adjust_bonf')
res=merge(res,da_p_c,all.x=TRUE)
res=res[order(-res$Count_down),]

res_1=res[res$Count_down>=12,]
#nrow(res_1)=600

res_2=res_1
for (i in (24+3):(24+2+24)){
	res_2=res_2[res_2[,i]<0 | is.na(res_2[,i]),]
}
#nrow(res_2)=282

res_3=res_2[res_2$CNV_corr_p_adjust_bonf<0.05 | is.na(res_2$CNV_corr_p_adjust_bonf),]
nrow(res_3)
#nrow(res_3)=266
res_3$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_3$peak)

to_plot=table(res_3$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('DOWN_Tumor_specific_DA_peaks.After_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='DOWN peaks after CNV correction',ylim=c(0,150))
dev.off()

res_2$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_2$peak)
to_plot=table(res_2$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('DOWN_Tumor_specific_DA_peaks.BeforeCNV_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='DOWN peaks before CNV correction',ylim=c(0,150))
dev.off()

peaks_1=StringToGRanges(res_3$peak, sep = c("-", "-"))
 ###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks=res_3
peaks$Gene=anno$SYMBOL
peaks$Type=anno$annotation
peaks$geneId=anno$geneId
peaks$peak_distanceToTSS=anno$distanceToTSS

#nrow(peaks[peaks$Type=='Promoter',])=2,907 2,907/3,829=76%

table(peaks$chr_peak[peaks$Type=='Promoter'])

write.table(peaks, "out/DOWN_Tumor_vsPT.Filtered.CNV_corrected.Annotated.20210811",sep='\t',row.names=F,quote=F)



#############################################################################################
#############################################################################################
#############################################################################################

#Karyotype plot:

regions=StringToGRanges(rownames(up_1[up_1$Count_down>=5,]), sep = c("-", "-"))
pdf("DOWN_in_atLeast_5_Tumors_vs_PT.pdf",useDingbats=F, height=7, width=7)
kp <- plotKaryotype(genome="hg38")
kpPlotRegions(kp, regions, col="blue")	
dev.off()