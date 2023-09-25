library(ggrepel)


####There are -log10(p-values) only in the Bulk ATAC-paper, and we calculate Z-score on their lo2(-log10(p_values))
setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/11.Compare_with_Bulk/')

me_bulk=read.table('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Archieved/Motif_analysis/Bulk_ATAC_ME.txt',header=TRUE,sep="\t")
me_bulk=me_bulk[!duplicated(me_bulk$TF_Name),]
rownames(me_bulk)=me_bulk$TF_Name
me_bulk=me_bulk %>% select ('TF_Name','Cluster_1')

#score=read.table('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/FindMarkers/2020-09-14/Motif_score_perMajorCell_group.tsv',header=TRUE,sep="\t")
score=read.table('Score_difference.cellGroups.20210617.tsv',header=TRUE,sep="\t")
score=score[,1:3]
colnames(score)=c('cell_type','TF_Name','mean_score')
score=score[score$cell_type=="Tumor",]

#score=read.table('Desktop/Projects/ATAC/CCRCC/Analysis/FindMarkers/2020-09-12/Manual_reannotation/out/Tumor cells_chromvar_MergedObj_motifs.tsv',header=TRUE,sep="\t")
#colnames(score)[7]='TF_Name'



score_1=merge(score, me_bulk)

score_1$mean_score_s=scale(score_1$mean_score)
score_1$Cluster_1_s=scale(log2(score_1$Cluster_1+1))


options(ggrepel.max.overlaps = Inf)
f_ch_cut=0.7
#degs=tab
score_1_sel=score_1[score_1$mean_score_s>f_ch_cut & score_1$Cluster_1_s>f_ch_cut,]

score_2_sel=score_1[score_1$mean_score_s>f_ch_cut & score_1$Cluster_1_s<f_ch_cut,]

p <- ggplot(score_1, aes(mean_score_s, Cluster_1_s)) + geom_point(alpha=0.5, colour="grey", data = score_1)

p <- p + geom_point(data = score_1, colour=ifelse(score_1$mean_score_s>f_ch_cut & score_1$Cluster_1_s>f_ch_cut,"blue","NA"), alpha=0.95)

p <- p + theme_minimal() #+geom_vline(xintercept=0, alpha=0.5)

p <- p + ylab("Motif Z-score bulk ATAC (Corces et al.)") 

p <- p + xlab("Motif Z-score snATAC-data") 

p <- p + geom_text_repel(data=score_2_sel,aes(y=Cluster_1_s,x= mean_score_s,label=TF_Name),alpha=1,colour = "black",segment.size=0.3)

p <- p + geom_text_repel(data=score_1_sel,aes(y=Cluster_1_s,x= mean_score_s,label=TF_Name),alpha=1,colour = "black",segment.size=0.3)

p <- p + geom_point(data = score_2_sel, colour=ifelse(score_2_sel$mean_score_s>f_ch_cut & score_2_sel$Cluster_1_s<f_ch_cut,"purple","NA"), alpha=0.95)

p <- p + xlim(-1,5)+ylim(-1,5)+geom_vline(xintercept=f_ch_cut, size=0.4,alpha=0.7,linetype=2)

p <- p + geom_hline(yintercept=f_ch_cut, size=0.4,alpha=0.7,linetype=2)



pdf(paste("ccRCC_snATAC_vsBulk_TopMotifs.2021-06-28.pdf",sep=""), width=7, height=7,useDingbats=FALSE)
print(p)
dev.off()
png(paste("Desktop/Projects/ATAC/CCRCC/Analysis/Compare_withBulk/2020-09-15/Kidney_bulk_snATAC_2020-09-16.png",sep=""), width=500, height=500)
p
dev.off()
