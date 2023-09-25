###Using mean_difference instead:
library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
library(ggrepel)


setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/17.Compare_withBulk/2021-10-08_ccRCCpeaks_atac/')

t1=read.table('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2021-04-09_ccRCC_peaks_ATAC/All_clusters_default.20210514/TFmotifView_enrichment_table_2021-05-14_19-56.txt',sep='\t',header=TRUE)
t2=read.table('Additional_motifs.20211008/TFmotifView_enrichment_table_2021-10-08_15-32.txt',sep='\t',header=T)
add_motifs=gsub('(.*)\\..*\\..*','\\1',t2$motif)
t2=t2[!(t2$motif %in% t1$motif),]
tab=rbind(t1,t2)
t3=read.table('../2021-10-08_mapPT_specificTFs/Additional_motifs.20211008/TFmotifView_enrichment_table_2021-10-08_18-10.txt',sep='\t',header=T)
t3=t3[!(t3$motif %in% tab$motif),]
tab=rbind(tab,t3)


degs=tab
degs$motif=gsub('(.*)\\..*\\..*','\\1',degs$motif)
degs=degs[!is.na(degs$global_fold_change),]
degs$FDR=p.adjust(degs$global_pvalue,method='fdr')
####Change FDR for two outliers (and malke gap in the figure)
degs$NegLog10FDR=-log10(degs$FDR)
degs$NegLog10FDR[1:2]=150

degs_sel=degs[degs$FDR<0.05,]
degs_sel=degs_sel[order(-degs_sel$global_fold_change),]
motifs_val=add_motifs[add_motifs %in% degs_sel$motif]
motifs_Notval=add_motifs[!(add_motifs %in% degs_sel$motif)]
degs_sel_1=degs[degs$FDR<0.05 & log2(degs$global_fold_change)>0.6,]


##########################
###Selected_tfs
##########################
options(ggrepel.max.overlaps = Inf)
genes=degs_sel$motif
genes_bulk=c(degs_sel$motif,add_motifs)
#genes=genes_bulk[genes_bulk %in% genes_sn]
genes=genes_bulk
genes_1=add_motifs
genes_2=degs_sel_1$motif[!(degs_sel_1$motif %in% add_motifs)]

cols=c(RColorBrewer::brewer.pal(n = 12, name = "Paired"),'black')
names(cols)=unique(as.character(tab$cell_t1))

p <- ggplot(degs, aes(log2(global_fold_change), NegLog10FDR)) + geom_point_rast(alpha=0.4, colour="grey", data = degs)

p <- p + geom_point(data=degs[degs$motif %in% c(motifs_val),],color='#800026', alpha=0.95,size=4)

p <- p + geom_point(data=degs[degs$motif %in% c(motifs_Notval),],color='blue', alpha=0.95,size=2)
 
p <- p + theme_minimal() +geom_hline(yintercept=-log10(0.05), alpha=0.5)+geom_vline(xintercept=0, alpha=0.5)

p <- p + ylab("-log10(FDR)") 

p <- p + xlab("Buk ATAC-seq Log2 Fold Change") 

p <- p + geom_text_repel(data=degs[degs$motif %in% c(add_motifs),],aes(y=NegLog10FDR,x= log2(global_fold_change),label=(motif)),alpha=1,colour = "black",segment.size=0.1,size=3.5)

p <- p + theme(strip.text.x = element_text(size = 12),  strip.text.y = element_text(size = 12))+
theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#+ylim(-50,370)+xlim(-6.5,6.5)


pdf(paste("Bulk_KIRC_specific_TFs.v3.20211008.pdf",sep=""), width=6, height=6.5,useDingbats=FALSE)
print(p)
dev.off()	


####NOW do for our data:
tab=read.table('Score_difference.Tumor_Normal_comparison.20210427.tsv',sep='\t',header=TRUE)
degs=tab
degs$diff=degs$mean_score2-degs$mean_score1
colnames(degs)[c(2,6)]=c('motif','global_fold_change')

degs_sel=degs[degs$FDR<0.05 ,]
degs_sel=degs_sel[order(-degs_sel$global_fold_change),]

genes=c(degs_sel$motif[28:48])
degs_bulk=degs_sel$motif[degs_sel$motif %in% genes_bulk & degs_sel$global_fold_change>0]
genes=degs_bulk
genes_sn=genes

cols=c(RColorBrewer::brewer.pal(n = 12, name = "Paired"),'black')
names(cols)=unique(as.character(tab$cell_t1))

p <- ggplot(degs, aes(global_fold_change, -log10(FDR))) + geom_point_rast(alpha=0.4, colour="grey", data = degs)

p <- p + geom_point(data=degs[degs$motif %in% c(genes),],color='blue', alpha=0.95,size=2)

p <- p +  theme_minimal() +geom_hline(yintercept=-log10(0.05), alpha=0.5)+geom_vline(xintercept=0, alpha=0.5)

p <- p + ylab("-log10(FDR)") 

p <- p + xlab("snATAC-seq") 

p <- p + geom_text_repel(data=degs[degs$motif %in% c(genes),],aes(y=-log10(FDR),x= global_fold_change,label=(motif)),alpha=1,colour = "black",segment.size=0.1)

p <- p + theme(strip.text.x = element_text(size = 12),  strip.text.y = element_text(size = 12))

p <- p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())

p <- p + labs(title="ccRCC-specific TF motifs")

pdf(paste("ccRCC_specific_TFs.v2.pdf",sep=""), width=9, height=6.5,useDingbats=FALSE)
print(p)
dev.off()	
