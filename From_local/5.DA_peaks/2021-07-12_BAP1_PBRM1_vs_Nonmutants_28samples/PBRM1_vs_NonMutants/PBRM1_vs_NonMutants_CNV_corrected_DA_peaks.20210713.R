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

date='20210713'
setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/5.DA_peaks/2021-07-12_BAP1_PBRM1_vs_Nonmutants_28samples/PBRM1_vs_NonMutants/')

da_peaks=read_delim('out/da_peaks_PBRM1mut_vs_NonMutants.min.pct0.1.20210713.tsv',delim='\t')
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


da_peaks=read_delim('out/da_peaks_PBRM1mut_vs_NonMutants.min.pct0.1.20210713.tsv',delim='\t')
da_peaks=as.data.frame(da_peaks)
da_peaks=da_peaks %>% dplyr::select ('peak','Sample','avg_log2FC')
up=dcast(da_peaks,peak~Sample,value.var='avg_log2FC')
colnames(up)[2:(ncol(up))]=paste('NonSignif_avg_lnFC',colnames(up)[2:(ncol(up))],sep='_')

res=merge(up_1,up,all.x=TRUE)


###cnv_corrected
da_p_c=read_delim('out/DA_peaks_PBRM1mutants_vs_NonMutants_correctedbyCNV.20210713.tsv',delim='\t')
da_p_c=as.data.frame(da_p_c)
colnames(da_p_c)[4]=paste('CNV_corr_',colnames(da_p_c)[4],sep='')
da_p_c=da_p_c %>% dplyr::select ('peak','CNV_corr_p_adjust_bonf')
res=merge(res,da_p_c,all.x=TRUE)
res=res[order(-res$Count_up),]

res_1=res[res$Count_up>=5,]
nrow(res_1) #2,042

res_2=res_1
for (i in 12:20){
	res_2=res_2[res_2[,i]>0 | is.na(res_2[,i]),]
}
nrow(res_2) #584

res_3=res_2[res_2$CNV_corr_p_adjust_bonf<0.05 | is.na(res_2$CNV_corr_p_adjust_bonf),]
nrow(res_3)
nrow(res_3) #561
res_3$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_3$peak)

to_plot=table(res_3$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('UP_PBRM1_specific_DA_peaks.After_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='UP peaks after CNV correction',ylim=c(0,100))
dev.off()

res_2$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_2$peak)
to_plot=table(res_2$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('UP_PBRM1_specific_DA_peaks.BeforeCNV_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='UP peaks before CNV correction',ylim=c(0,100))
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

#nrow(peaks[peaks$Type=='Promoter',])= .../561=...%

table(peaks$chr_peak[peaks$Type=='Promoter'])

write.table(peaks, "out/UP_PBRM1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210713.tsv",sep='\t',row.names=F,quote=F)


---------------------------------
###############
#NOW DOWN PEAKS
###############
---------------------------------
da_peaks=read_delim('out/da_peaks_PBRM1mut_vs_NonMutants.min.pct0.1.20210713.tsv',delim='\t')
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


da_peaks=read_delim('out/da_peaks_PBRM1mut_vs_NonMutants.min.pct0.1.20210713.tsv',delim='\t')
da_peaks=as.data.frame(da_peaks)
da_peaks=da_peaks %>% dplyr::select ('peak','Sample','avg_log2FC')
up=dcast(da_peaks,peak~Sample,value.var='avg_log2FC')
colnames(up)[2:(ncol(up))]=paste('NonSignif_avg_lnFC',colnames(up)[2:(ncol(up))],sep='_')

res=merge(up_1,up,all.x=TRUE)


###cnv_corrected
da_p_c=read_delim('out/DA_peaks_PBRM1mutants_vs_NonMutants_correctedbyCNV.20210713.tsv',delim='\t')
da_p_c=as.data.frame(da_p_c)
colnames(da_p_c)[4]=paste('CNV_corr_',colnames(da_p_c)[4],sep='')
da_p_c=da_p_c %>% dplyr::select ('peak','CNV_corr_p_adjust_bonf')
res=merge(res,da_p_c,all.x=TRUE)
res=res[order(-res$Count_down),]

res_1=res[res$Count_down>=5,]
#nrow(res_1)=444

res_2=res_1
for (i in 12:20){
	res_2=res_2[res_2[,i]<0 | is.na(res_2[,i]),]
}
#nrow(res_2)=87

res_3=res_2[res_2$CNV_corr_p_adjust_bonf<0.05 | is.na(res_2$CNV_corr_p_adjust_bonf),]
nrow(res_3)
#nrow(res_3)=85
res_3$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_3$peak)

to_plot=table(res_3$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('DOWN_PBRM1_specific_DA_peaks.After_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='UP peaks after CNV correction',ylim=c(0,60))
dev.off()

res_2$chr_peak=gsub('chr(.*)-.*-.*','\\1',res_2$peak)
to_plot=table(res_2$chr_peak)
to_plot=to_plot[order(names(to_plot))]
to_plot=to_plot[order(factor(names(to_plot),levels=c(1:22,"X")))]
pdf(paste('DOWN_PBRM1_specific_DA_peaks.BeforeCNV_correction.',date,'.pdf',sep=''),width=8,height=7,useDingbats=F)
barplot(to_plot,ylab='Number of peaks',main='UP peaks before CNV correction',ylim=c(0,60))
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

#nrow(peaks[peaks$Type=='Promoter',])=52 52/85=62%

table(peaks$chr_peak[peaks$Type=='Promoter'])

write.table(peaks, "out/DOWN_PBRM1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210713.tsv",sep='\t',row.names=F,quote=F)


----------------------------
#Compare with the BAP1-specific
up_pbrm1=read_delim('out/UP_PBRM1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210623.tsv',delim='\t')
up_pbrm1=as.data.frame(up_pbrm1)

down_bap1=read_delim('../BAP1_vs_NonMutants/out/DOWN_BAP1_vsNonMutants.Filtered.CNV_corrected.Annotated.20210623.tsv',delim='\t')
down_bap1=as.data.frame(down_bap1)



-----------------------------

m1=cbind(c(134,647),c(2970,3890))

#it seems consistent with the results in the paper: https://www.nature.com/articles/s41467-018-08255-x (fig.2F)

> 598+651
[1] 1249
> 317
[1] 317
> 317/1249
[1] 0.253803
> 647/3890
[1] 0.1663239
 
