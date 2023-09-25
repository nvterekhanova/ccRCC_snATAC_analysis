###soource script: Desktop/Projects/ATAC/CCRCC/Analysis/Markers.20201205/Cell_groups_5C/20201206/Cell_group_markers_5C.20201206.R
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


theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())

setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/9.Motifs_per_cellGroup/20210809_Cell_group_motifs_5B/')

cell_types=c('Immune','Stroma','Tumor','PT')
score=read.table('Score_difference.cellGroups.20210808.tsv',header=TRUE,sep="\t")
score=score[score$cell_t1 %in% cell_types,]
score=score[,1:3]
colnames(score)=c('cell_type','TF_Name','mean_score')

########################################
####Trying to make only unique (our latest versions):
###############################


#Plotting all scores in a heatmap
score_diff=read.table('Score_difference.cellGroups.20210808.tsv',sep='\t',header=TRUE)
score_diff=score_diff[score_diff$cell_t1 %in% cell_types,]
score_diff=score_diff[order(-score_diff$diff),]
colnames(score_diff)[1]='cell_t2'

score_stat=score_diff[score_diff$mean_score1>0,]
score_stat=dcast(score_stat,TF_Name~cell_t2,value.var='FDR')
score_stat$na_count <- apply(score_stat[,-1], 1, function(x) sum(is.na(x)))
score_stat_Epit=score_stat[score_stat$Tumor<0.05 & score_stat$PT<0.05,]
#score_stat=score_stat[score_stat$na_count==4,]
score_stat_Epit=score_stat_Epit[!is.na(score_stat_Epit$TF_Name),]

score_diff=score_diff[order(-score_diff$mean_score1),]
score_diff=score_diff[!duplicated(score_diff$TF_Name),]
all_top=NULL
for (i in unique(as.character(score_diff$cell_t2))){
	if(i %in% c('Tumor','PT')){
	top=score_diff[score_diff$cell_t2==i & score_diff$FDR<0.05 & score_diff$mean_score1>0,][1:20,]
	}else if(i=='Stroma'){
	top=score_diff[score_diff$cell_t2==i & score_diff$FDR<0.05 & score_diff$mean_score1>0,]
	top=top[!grepl('FOS|JUN|BACH|JDP|NFE|BATF',top$TF_Name),][1:10,]
	}else{
	top=score_diff[score_diff$cell_t2==i & score_diff$FDR<0.05 & score_diff$mean_score1>0,][1:10,]
	}
	all_top=rbind(all_top,top)
}
all_top=rbind(all_top,score_diff[score_diff$cell_t2=="Tumor cells" & score_diff$TF_Name=="ARNT::HIF1A",])



markers=score[score$TF_Name %in% all_top$TF_Name,]
enr=all_top[,1:2]
colnames(enr)[1]="Cell_type_enriched"
markers=merge(markers,enr,all.x=TRUE)



cols <- brewer.pal(9, "YlOrRd")
getPalette= colorRampPalette(cols)

all_merged=markers
all_merged$Cell_type_enriched=as.character(unlist(all_merged$Cell_type_enriched))
all_merged$cell_type=as.character(unlist(all_merged$cell_type))


all_merged$cell_type=factor(all_merged$cell_type,levels=rev(c('Immune','Stroma','PT','Tumor')))
all_merged=all_merged[!is.na(all_merged$cell_type),]
all_merged$Cell_type_enriched=factor(all_merged$Cell_type_enriched,levels=rev(c('Immune','Stroma','PT','Tumor')))


####Heatmap-version:
all_merged$z_mean_score=scale(all_merged$mean_score)
#cols <- c("#377EB8","white","#E41A1C")
cols <- c("#377EB8","#377EB8","white","#E41A1C","#E41A1C")
getPalette= colorRampPalette(cols)

print(min(all_merged$mean_score))
print(max(all_merged$mean_score))

all_merged$mean_score=ifelse(all_merged$mean_score>5,5,all_merged$mean_score)
all_merged$mean_score=ifelse(-all_merged$mean_score>5,-5,all_merged$mean_score)

p <-  ggplot()

p <- p + facet_grid(.~Cell_type_enriched,drop=T,scales = "free", space = "free")

p <- p + geom_tile(data=all_merged,aes(x=TF_Name, y=cell_type, fill= mean_score), linetype="blank",width=0.9, height=0.9)

p <- p + scale_fill_gradientn(name= "Motif Score", colours=getPalette(100), na.value="grey", limit=c(-5,5))

p <- p  + theme_bw() + theme_nogrid()

p <- p + theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.2,hjust=0.95), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())

p <- p + theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))

p <- p + theme(strip.text.x=element_text(size=12)
)


pdf(paste("Specific_motifs_cellGroups.v1.2021-08-09.pdf",sep=""),width=12, height=3.1,useDingbats=FALSE)
p
dev.off()