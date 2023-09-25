####the latest version of the plot (include more categories)
####Making summary of peaks:
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Signac)


all_peaks=NULL
samples=c('')
for (sample in samples)	{
	tab=read.table(paste("/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/QC_summary/peaks_summary.20210829/recentered_peaks/recentered_final.filtered",sample,".tsv",sep=""),header=TRUE,sep="\t")
	peaks_1=StringToGRanges(tab$new_peak, sep = c("-", "-"))
	peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
	anno=as.data.frame(peakAnno@anno)
    	anno$annotation=gsub('(Exon).*','\\1',anno$annotation)
    	anno$annotation=gsub('(Intron).*','\\1',anno$annotation)
	anno$annotation=gsub('Downstream.*','Downstream (<3kb)',anno$annotation)
	tab$peak_type=anno$annotation
	tab$Count=1
	tab_1=aggregate(tab$Count, by=list(tab$peak_type),FUN='sum')
	colnames(tab_1)=c('Type','Count')
	tab_1$Sample=sample
	all_peaks=rbind(all_peaks,tab_1)
	print(sample)
}

sample_data=read.table("/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/DATA/meta_data.20200505.v1.tsv",header=TRUE,sep="\t")
s_data=sample_data %>% dplyr::select ("Aliquot.snRNA","Aliquot.snRNA.WU")
colnames(s_data)[1]="Sample"
all_s1=merge(all_peaks, s_data,all.x=TRUE)

all_s1$Sample_type=ifelse(all_s1$Sample %in% c("CPT0075170013","CPT0000890002","CPT0014470002", "CPT0001270002"),"Normal","Tumor")
all_s1$Sample_type=factor(all_s1$Sample_type,levels=c('Tumor','Normal')) 

all_s1$Type=factor(all_s1$Type,levels=c("Promoter","5' UTR","Exon","Intron","3' UTR","Downstream (<3kb)","Distal Intergenic")) 

###for summary numbers in the text:
st=aggregate(all_s1$Count, by=list(all_s1$Sample),FUN='sum')


all_s1=all_s1 %>% mutate(Aliquot.snRNA.WU =  factor(Aliquot.snRNA.WU, levels = order)) %>% arrange(Aliquot.snRNA.WU)


write.table(all_s1,"/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/QC_summary/peaks_summary.20210829/Peaks_annotation_acrossSamples_forSummaryFig.20210902.tsv",sep='\t',quote=F,row.names=F)
#cols=c('promoter'='#E41A1C','intergenic'='#377EB8','distal'='#4DAF4A')
cols=brewer.pal(n = 7, name = "Set1")

p <- ggplot(all_s1, aes(x = Type, y = Count)) + 

p  <- p + geom_bar(stat="identity",aes(fill = Type))

p <- p + facet_grid(.~Sample_type+Aliquot.snRNA.WU,scales = "free", space = "free") 

p <- p + theme(axis.text.x = element_text(colour="black", size=10, angle=45, vjust = 1,hjust=1),
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank()) + labs(x="")

p <- p + scale_fill_manual(values=cols) + theme_minimal()

p <- p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),
axis.ticks = element_blank(),strip.text.x=element_text(size=12,angle=90,vjust=0.2,hjust=0.95))

p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- p + theme(legend.position='bottom', legend.text=element_text(size=12))

pdf(paste("/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/QC_summary/peaks_summary.20210829/peaks_summary_2021-09-02.pdf",sep=""), width=14, height=5,useDingbats=FALSE)
p
dev.off()