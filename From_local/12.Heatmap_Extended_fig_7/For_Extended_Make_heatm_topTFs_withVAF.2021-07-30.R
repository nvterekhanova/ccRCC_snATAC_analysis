library(Signac)
library(Repitools)
library(ChIPseeker)
library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(reshape)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)

ecdf_fun = function(x,perc) ecdf(x)(perc)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())


setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/12.Heatmap_Extended_fig_7/')

###Plot all scores:
top_tfs=c('HIF1A','NFKB1','NFKB2')
tfs_sel=top_tfs
expr=read.table('Motif_score_perCell_group.20210715.tsv',sep='\t',header=TRUE)
expr$cell_type=gsub('(.*)_.*','\\1',expr$cell_type)

all_2=dcast(expr,cell_type~TF_Name,value.var="mean_score")
all_2=all_2[,colnames(all_2) %in% c('cell_type',top_tfs)]
expr=all_2
rownames(expr)=expr[,1]
expr=expr[,-1]
expr=t(expr)
expr=as.data.frame(expr)

tab=expr
selected_ids=gsub('(.*)-(.*)-(.*)[0-9]','\\1-\\2-\\3',colnames(tab))
g=c('BAP1','PBRM1','SETD2','KDM5C')
bulk_rna=read.table("/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Archieved/Motif_analysis/2020-06-07/bulk_DEGs/CCRCC_RNA_rpkm_tumor_normal_ids.tsv",header=TRUE,sep="\t")
bulk_rna=bulk_rna[!duplicated(bulk_rna$geneID),]
row.names(bulk_rna)=bulk_rna[,1]
bulk_rna=bulk_rna[,-1]
bulk_rna=log2(bulk_rna+1)
bulk_rna_scaled=t(bulk_rna)
bulk_rna_scaled=scale(bulk_rna_scaled)
bulk_rna_scaled=t(bulk_rna_scaled)
bulk_rna=bulk_rna_scaled
colnames(bulk_rna)=gsub('\\.','-',colnames(bulk_rna))
bulk_rna=as.data.frame(bulk_rna)
colnames(bulk_rna)=gsub('_Normal','-N',colnames(bulk_rna))
colnames(bulk_rna)=gsub('_Tumor','-T',colnames(bulk_rna))
bulk_rna=bulk_rna[,selected_ids]
col_rna=bulk_rna[rownames(bulk_rna) %in% g,]
colnames(col_rna)=colnames(tab)



protein=read.table('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/FindMarkers/2020-09-29/proteome_per_gene_CCRCC_Tumor_Normal.tsv',sep='\t',header=TRUE)
samples=read.table('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/FindMarkers/2020-09-29/Protein_sample_data.tsv',sep='\t',header=TRUE)
samples$Type=ifelse(as.character(unlist(samples$sample_type))=="Primary Tumor","Tumor","Normal")
samples$ID=paste(samples$case_id,"_",samples$Type,sep='')
samples=samples[colnames(protein),]
colnames(protein)=samples$ID
colnames(protein)=gsub('_Normal','-N',colnames(protein))
colnames(protein)=gsub('_Tumor','-T',colnames(protein))
protein=protein[,selected_ids]
protein=protein[rownames(protein) %in% g,]
colnames(protein)=colnames(tab)
prot_scale=t(protein)
prot_scale=scale(prot_scale)
protein=t(prot_scale)
col_protein=as.data.frame(protein)


col_rna=t(col_rna)
col_protein=t(col_protein)
col_rna=as.data.frame(col_rna)
col_protein=as.data.frame(col_protein)
colnames(col_rna)=paste('RNA','_',colnames(col_rna)[1:4],sep='')


meta=read.table('Sample_categories.20210503.txt',sep='\t',header=T)
bap1_s=meta$Aliquot.WU[meta$Category=='BAP1-mutant']
pbrm1_s=meta$Aliquot.WU[meta$Category=='PBRM1-mutant']
both_s=meta$Aliquot.WU[meta$Category=='Both mutated']
bap1_s=bap1_s[bap1_s!='']
nat=c('')
nonm_s=c(meta$Aliquot.WU[meta$Category=='Non-mutant'],'')

all=cbind(col_rna,col_protein)
all$Sample=row.names(all)
all$Type=gsub('.*_(.*)','\\1',all$Sample)
all$Type=ifelse(all$Type=="T2","Tumor",all$Type)
all$BAP1_mutant=ifelse(all$Sample %in% c(bap1_s,both_s),1,0)
all$PBRM1_mutant=ifelse(all$Sample %in% c(pbrm1_s,both_s),1,0)
all$case=gsub('(.*-.*)-.*','\\1',all$Sample)

vaf=read.table('PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv',sep='\t',header=T)
vaf=vaf[,colnames(vaf) %in% c('Case',g)]
vaf1=as.data.frame(lapply(vaf, function(x) as.numeric(gsub('.*\\((.*)\\)','\\1',x))))
vaf1$Case=vaf$Case
rownames(vaf1)=vaf1$Case
vaf1=vaf1[all$case,]
vaf1$Sample=all$Sample
vaf1$BAP1[vaf1$Sample=='']<-NA


all$BAP1_mutant=vaf1$BAP1
all$PBRM1_mutant=vaf1$PBRM1
all$KDM5C_mutant=vaf1$KDM5C
all$SETD2_mutant=vaf1$SETD2

color_rna=colorRamp2(c(-2, 0, 2), c("#0571b0", "white", "#ca0020"))
all$Type=ifelse(all$Sample %in% nat,'NAT','Tumor')
split=factor(all$Type,levels=c('NAT','Tumor'))

column_ha = HeatmapAnnotation(BAP1_Mutant=all$BAP1_mutant,PBRM1_mutant=all$PBRM1_mutant,KDM5C_mutant=all$KDM5C_mutant,SETD2_mutant=all$SETD2_mutant,RNA_BAP1=all$RNA_BAP1, RNA_PBRM1=all$RNA_PBRM1,RNA_KDM5C=all$RNA_KDM5C,RNA_SETD2=all$RNA_SETD2, Protein_BAP1=all$BAP1, Protein_PBRM1=all$PBRM1,Protein_KDM5C=all$KDM5C,Protein_SETD2=all$SETD2,col=list(BAP1_Mutant=colorRamp2(c(0, 0.5), c("white smoke", '#984EA3')),PBRM1_mutant=colorRamp2(c(0, 0.5), c("white smoke", '#FF7F00')), SETD2_mutant=colorRamp2(c(0, 0.5), c("white smoke", '#377eb8')),KDM5C_mutant=colorRamp2(c(0, 0.5), c("white smoke", '#a65628')),RNA_BAP1=color_rna,RNA_PBRM1=color_rna,RNA_KDM5C=color_rna,RNA_SETD2=color_rna,Protein_BAP1=color_rna,Protein_PBRM1=color_rna,Protein_KDM5C=color_rna,Protein_SETD2=colorRamp2(c(-2, 0, 2), c("#0571b0", "white", "#ca0020"))),gp = gpar(col = "black"),show_legend = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE))

tab=tab[top_tfs,]

col_order=c(nat,nonm_s,pbrm1_s,both_s,bap1_s)
h=Heatmap(tab,top_annotation = column_ha,column_split=split,show_parent_dend_line = FALSE,cluster_rows = FALSE,cluster_columns=FALSE,column_order=col_order)

pdf(paste("For_the_ExtendedFig.9a.VAFadded.20210729.pdf",sep=""), width=10, height=5,useDingbats=FALSE)
print(h)
dev.off()
	