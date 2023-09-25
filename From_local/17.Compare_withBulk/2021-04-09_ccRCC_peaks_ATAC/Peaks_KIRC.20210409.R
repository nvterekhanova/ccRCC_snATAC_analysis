####2020-06-23:
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Signac)
library(Seurat)
libraary(ComplexHeatmap)
library(circlize)
library(reshape)
library(reshape2)
library(tidyverse)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

###https://benbermanlab.com/assets/code/Workshop%20for%20ATAC-seq%20analysis.html
setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-04-09_ccRCC_peaks_ATAC/')

###these are just peak-coordinates:
peaks=read_delim('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-03-12_KIRC_bulk_ATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/KIRC_peakCalls.txt',delim='\t')
peaks=as.data.frame(peaks)

####Try log2 matrices for all cancer types
peaks=read_delim('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-04-09_ccRCC_peaks_ATAC/TCGA-ATAC_PanCan_Log2Norm_Counts.txt',delim='\t')
peaks=as.data.frame(peaks)
peaks_1=peaks

p=readRDS('/Users/nadezhdaterekhanova/Downloads/TCGA-ATAC_PanCan_Log2Norm_Counts.rds')
p=as.data.frame(p)
peaks=p

samples.ids=read_delim('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-03-12_KIRC_bulk_ATAC/TCGA_identifier_mapping.txt',delim='\t')
samples.ids=as.data.frame(samples.ids)

#colnames(peaks)[-c(1:5)] <- samples.ids$Case_ID[match(gsub("_","-",colnames(peaks)[-c(1:5)]),samples.ids$bam_prefix)]
#colnames(peaks)[-c(1:5)]=gsub('(.*-.*-.*-.*)-.*-.*-.*','\\1',colnames(peaks)[-c(1:5)])
samples.ids$Cancer=gsub('(.*)-.*-.*-.*-.*-.*-.*-.*-.*-.*-.*-.*','\\1',samples.ids$bam_prefix)
samples.ids$BAM_p=gsub('-','_',samples.ids$bam_prefix)


kirc=samples.ids$BAM_p[samples.ids$Cancer=='KIRC']
kirp=samples.ids$BAM_p[samples.ids$Cancer=='KIRP']
other=samples.ids$BAM_p[!(samples.ids$Cancer %in% c('KIRC'))]

all_st=NULL
for (n in 0:563){
all_st_1=NULL
for(i in ((n*1000+1):((n+1)*1000))){
	t1_p=as.numeric(peaks[i,colnames(peaks) %in% kirc])
	t2_p=as.numeric(peaks[i,colnames(peaks) %in% c(other)])
	f_ch=mean(t1_p)-mean(t2_p)
	test=t.test(t1_p,t2_p,alternative=c('two.sided'))
#	test=t.test(t1_p,t2_p,alternative=c('greater'))
	st= cbind(peaks[i,1:5],test$p.value,f_ch)
	all_st_1=rbind(all_st_1,st)
	print(i)
	}
	all_st=rbind(all_st,all_st_1)
}
all_st=rbind(all_st,all_st_1)
all_st=as.data.frame(all_st)
colnames(all_st)[6:7]=c('p_value','LnFch')
all_st$FDR=p.adjust(all_st$p_value,method='fdr')
all_st=all_st[order(all_st$FDR),]
write.table(all_st,'Bulk_ATAC_peaks/KIRC_vsOthers.v2.tsv',sep='\t',row.names=FALSE,quote=FALSE)

###################################################
###Now look for peaks in comparison KIRC vs KIRP###
###################################################
all_st=NULL
for (n in 0:563){
all_st_1=NULL
for(i in ((n*1000+1):((n+1)*1000))){
	t1_p=as.numeric(peaks[i,colnames(peaks) %in% kirc])
	t2_p=as.numeric(peaks[i,colnames(peaks) %in% kirp])
	f_ch=mean(t1_p)-mean(t2_p)
	test=t.test(t1_p,t2_p,alternative=c('two.sided'))
#	test=t.test(t1_p,t2_p,alternative=c('greater'))
	st= cbind(peaks[i,1:5],test$p.value,f_ch)
	all_st_1=rbind(all_st_1,st)
	print(i)
	}
	all_st=rbind(all_st,all_st_1)
}
all_st=rbind(all_st,all_st_1)
all_st=as.data.frame(all_st)
colnames(all_st)[6:7]=c('p_value','LnFch')
all_st$FDR=p.adjust(all_st$p_value,method='fdr')
all_st=all_st[order(all_st$FDR),]
write.table(all_st,'Bulk_ATAC_peaks/KIRC_vsKIRP.v2.tsv',sep='\t',row.names=FALSE,quote=FALSE)

##################################################
##################################################
##################################################
#Then compare, find overlap btween two sets of peaks:
setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-04-09_ccRCC_peaks_ATAC/')
kirc_other=read.table('Bulk_ATAC_peaks/KIRC_vsOthers.v2.tsv',sep='\t',header=TRUE)
kirc_kirp=read.table('Bulk_ATAC_peaks/KIRC_vsKIRP.v2.tsv',sep='\t',header=TRUE)
colnames(kirc_other)[6:8]=paste('KIRC_OTHER',colnames(kirc_other)[6:8],sep='_')
colnames(kirc_kirp)[6:8]=paste('KIRC_KIRP',colnames(kirc_kirp)[6:8],sep='_')
both=merge(kirc_other,kirc_kirp)
both=both[both$KIRC_OTHER_FDR<0.05 & both$KIRC_KIRP_FDR<0.05,]
both1=both[both$KIRC_OTHER_LnFch>0 & both$KIRC_KIRP_LnFch>0,]
#both=both[(both$KIRC_OTHER_LnFch>0.6 & both$KIRC_KIRP_LnFch>0.6),]
write.table(both1,'Bulk_ATAC_peaks/Overlap_KIRC_specific.v4.tsv',sep='\t',row.names=FALSE,quote=FALSE)
bed=both1[,c('seqnames','start','end')]
colnames(bed)=c('chromosome','starting position','ending position')
write.table(bed,'Bulk_ATAC_peaks/Overlap_KIRC_specific.v4.bed',sep='\t',row.names=FALSE,quote=FALSE)





all_st=NULL
for(i in 1:nrow(peaks)){
	t1_p=as.numeric(peaks[i,colnames(peaks) %in% pbrm1])
	t2_p=as.numeric(peaks[i,colnames(peaks) %in% non_m])
	test=t.test(t1_p,t2_p,alternative=c('greater'))
	st= cbind(peaks[i,1:5],test$p.value)
	all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
colnames(all_st)[6]='p_value'
all_st$FDR=p.adjust(all_st$p_value,method='fdr')
all_st=all_st[order(all_st$FDR),]
write.table(all_st,'peaks_BAP1_PBRM1/PBRM1_up.tsv',sep='\t',row.names=FALSE,quote=FALSE)

all_st_s=all_st[all_st$FDR<0.05,]
write.table(all_st_s,'peaks_BAP1_PBRM1/PBRM1_up.FDR_0.05.tsv',sep='\t',row.names=FALSE,quote=FALSE)

all_st=read.table('peaks_BAP1_PBRM1/BAP1_up.tsv',sep='\t',header=TRUE)
all_st_s=all_st[all_st$FDR<0.05,]
write.table(all_st_s,'peaks_BAP1_PBRM1/BAP1_up.FDR_0.05.tsv',sep='\t',row.names=FALSE,quote=FALSE)

all_st=NULL
for(i in 1:nrow(peaks)){
	t1_p=as.numeric(peaks[i,colnames(peaks) %in% non_m])
	t2_p=as.numeric(peaks[i,colnames(peaks) %in% pbrm1])
	test=t.test(t1_p,t2_p,alternative=c('greater'))
	st= cbind(peaks[i,1:5],test$p.value)
	all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
colnames(all_st)[6]='p_value'
all_st$FDR=p.adjust(all_st$p_value,method='fdr')
all_st=all_st[order(all_st$FDR),]
write.table(all_st,'peaks_BAP1_PBRM1/Non_m_vsPBRM1_up.tsv',sep='\t',row.names=FALSE,quote=FALSE)

all_st_s=all_st[all_st$FDR<0.05,]
write.table(all_st_s,'peaks_BAP1_PBRM1/Non_m_vsPBRM1_up.FDR_0.05.tsv',sep='\t',row.names=FALSE,quote=FALSE)



####Now do some testing:
###BAP1-mutants vs others:
####Try log2 matrices for all cancer types
peaks=read_delim('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-04-09_ccRCC_peaks_ATAC/TCGA-ATAC_PanCan_Log2Norm_Counts.txt',delim='\t')
peaks=as.data.frame(peaks)
peaks_1=peaks

p=readRDS('/Users/nadezhdaterekhanova/Downloads/TCGA-ATAC_PanCan_Log2Norm_Counts.rds')
p=as.data.frame(p)
peaks=p

samples.ids=read_delim('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-03-12_KIRC_bulk_ATAC/TCGA_identifier_mapping.txt',delim='\t')
samples.ids=as.data.frame(samples.ids)

#colnames(peaks)[-c(1:5)] <- samples.ids$Case_ID[match(gsub("_","-",colnames(peaks)[-c(1:5)]),samples.ids$bam_prefix)]
#colnames(peaks)[-c(1:5)]=gsub('(.*-.*-.*-.*)-.*-.*-.*','\\1',colnames(peaks)[-c(1:5)])
samples.ids$Cancer=gsub('(.*)-.*-.*-.*-.*-.*-.*-.*-.*-.*-.*-.*','\\1',samples.ids$bam_prefix)
samples.ids$BAM_p=gsub('-','_',samples.ids$bam_prefix)
samples.ids$Case_id_2=gsub('(.*-.*-.*-.*)-.*-.*-.*','\\1',samples.ids$Case_ID)


kirc=samples.ids$BAM_p[samples.ids$Cancer=='KIRC']
bap1=samples.ids$BAM_p[samples.ids$Case_id_2=='TCGA-B8-A54G-01A']
not_bap1=kirc[!(kirc %in% bap1)]
#kirp=samples.ids$BAM_p[samples.ids$Cancer=='KIRP']
#other=samples.ids$BAM_p[!(samples.ids$Cancer %in% c('KIRC'))]
#tab=peaks[,colnames(peaks) %in% kirc]

#bap1_s=c('TCGA-B8-A54G-01A')

all_st=NULL
for (n in 0:563){
all_st_1=NULL
for(i in ((n*1000+1):((n+1)*1000))){
	t1_p=as.numeric(peaks[i,colnames(peaks) %in% bap1])
	t2_p=as.numeric(peaks[i,colnames(peaks) %in% not_bap1])
	f_ch=mean(t1_p)-mean(t2_p)
	test=t.test(t1_p,t2_p,alternative=c('two.sided'))
#	test=t.test(t1_p,t2_p,alternative=c('greater'))
	st= cbind(peaks[i,1:5],test$p.value,f_ch)
	all_st_1=rbind(all_st_1,st)
	print(i)
	}
	all_st=rbind(all_st,all_st_1)
}
all_st=rbind(all_st,all_st_1)
all_st=as.data.frame(all_st)
colnames(all_st)[6:7]=c('p_value','LnFch')
all_st$FDR=p.adjust(all_st$p_value,method='fdr')
all_st=all_st[order(all_st$FDR),]
write.table(all_st,'BAP1_specific/BAP1_vsOtherKIRC.20210610.tsv',sep='\t',row.names=FALSE,quote=FALSE)


all_st_f=all_st[all_st$FDR<0.05,]
all_st_f$peak=paste(all_st_f$seqnames,all_st_f$start,all_st_f$end,sep='-')


peaks_1=StringToGRanges(all_st_f$peak, sep = c("-", "-"))
 ###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks_1=all_st_f
peaks_1$Gene=anno$SYMBOL
peaks_1$Type=anno$annotation
peaks_1$Cat=gsub('(.*)_.*','\\1',rownames(peaks_1))

down=peaks_1[peaks_1$LnFch<0,]
up=peaks_1[peaks_1$LnFch>0,]


