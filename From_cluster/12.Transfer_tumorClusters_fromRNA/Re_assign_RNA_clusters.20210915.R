library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(reshape2)
library(RColorBrewer)

library(EnsDb.Hsapiens.v86)
args = commandArgs(trailingOnly=TRUE)

samples=c('')

all_res_sel=NULL
all_tab=NULL
for (sample in samples){
print(sample)

ATAC=readRDS(paste("TumorOnly.RDS/",sample,"_TumorOnly.20210913.rds",sep=""))
clusters=read.table(paste('TumorOnly_RNAClusters_mapped/metadata/',sample,'_TumorOnly.snATAC.meta.data',sep=''),
sep='\t',header=T)

###Read RNA-object:
rna_paths=read.table('data_snRNA/Paths_TumorCellOnlyReclustered_SeuratObject.20210805.v1.tsv',
sep='\t',header=T)
rna_paths=rna_paths[rna_paths$Aliquot.snRNA %in% samples,]
path_rna=rna_paths$Path_katmai[rna_paths$Aliquot.snRNA==sample]

RNA=readRDS(path_rna)
case=rna_paths$Aliquot.snRNA.WU[rna_paths$Aliquot.snRNA==sample]

##############

x=as.data.frame(table(clusters %>% dplyr::select ('seurat_clusters','Cluster_ID')))
x1=aggregate(x$Freq, by=list(x$seurat_clusters),FUN='sum')
colnames(x1)=c('seurat_clusters','Sum')

res=merge(x,x1,all.x=T)
res$Fraction=res$Freq/res$Sum
res_sel=res[res$Fraction>0.5,]

orig_1=as.data.frame(ATAC$seurat_clusters)
colnames(orig_1)='seurat_clusters'
orig_1$Barcode=rownames(orig_1)
meta=orig_1
meta=merge(meta,res_sel,all.x=T)
rownames(meta)=meta$Barcode
meta=meta[rownames(orig_1),]
ATAC$RNA_Cluster_ID=meta$Cluster_ID

orig_1=as.data.frame(ATAC$RNA_Cluster_ID)
colnames(orig_1)='cell_group'
orig_1$sample_barcode=paste(sample,rownames(orig_1),sep='_')
orig_1$cell_group=paste(case,orig_1$cell_group,sep='_')
all_tab=rbind(all_tab,orig_1)



#####Labelling cell-types for the snRNA

meta=read.table('data_snRNA/Barcode_TumorPTLOH_ByCluster.20210903.v1.tsv',sep='\t',header=T)
meta$Sample=gsub('(.*)_.*','\\1',meta$sample_barcode)
meta_sample=meta[meta$Sample==sample,]
meta_sample$Cluster_group=gsub('(.*)_.*','\\1',meta_sample$cell_group)
meta_sample=meta_sample[meta_sample$Cluster_group==case,]
meta_sample$Cluster_ID=gsub('(.*)_(.*)','\\2',meta_sample$cell_group)
rownames(meta_sample)=gsub('(.*)_(.*)','\\2',meta_sample$sample_barcode)

orig_1=as.data.frame(RNA$orig.ident)
meta_sample=meta_sample[rownames(orig_1),]

RNA$RNA_Cluster_ID=meta_sample$Cluster_ID

plot1 <- DimPlot(
  object = RNA,
  group.by = 'RNA_Cluster_ID',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(case,' snRNA-seq',sep=""))

plot2 <- DimPlot(
  object = ATAC,
  group.by = 'RNA_Cluster_ID',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(case,' snATAC-seq',sep=""))

pdf(paste("TumorOnly.RDS.dimplots.RNAClustersReassigned/plots/",sample,"_",case,
"_snRNA_ClusterIDs_mapped_to_snATAC.pdf",sep=""),
height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()

res_sel$Case=case
res_sel$Sample=sample
all_res_sel=rbind(all_res_sel,res_sel)

}

write.table(all_tab,paste("TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/RNA_clusters_TumorOnly.snATAC.tsv",
sep=""),sep='\t',quote=FALSE)
write.table(all_res_sel,paste("TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/RNA_clusters_Fractions_",
"inTumorOnly_snATAC_Clusters.tsv",sep=""),sep='\t',quote=FALSE)


####2021-09-17, need to add PT-clusters here as well (selected RNA-based PT-clusters (0, 5, 8, 13, 4, 3, 11)

tumor_map=read.table('TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/RNA_clusters_TumorOnly.snATAC.tsv',
sep='\t',header=T)

nat=readRDS(paste('../10.Cell_of_origin/Map_snRNA_clusters.20210804/NAT_PT_LOH_cellsOnly.20210804.',
'v1.GeneActivity.rds.gz',sep=''))
nat_meta=nat@meta.data
cl=table(nat_meta$RNA_seurat_clusters[nat_meta$cell_type=='PT'])
pt_clusters_rna=names(cl[cl>200])

###we also won't use cluaters 12, 14, because they have mix of LOH and PT based on slides 39-40: https://docs.google.com/presentation/d/1dKCXGKxwPx4VaUfZL_5xGVlVyqO-xaNJ/edit#slide=id.p39


pt_clusters_rna=c(0,5,8,13,4,3,11)
nat_meta$ID=paste(nat_meta$cell_type,nat_meta$RNA_seurat_clusters,sep='_')
nat_meta$ID=ifelse(nat_meta$RNA_seurat_clusters %in% pt_clusters_rna & nat_meta$cell_type=='PT',
nat_meta$ID, 'NA')
nat_meta$sample_barcode=rownames(nat_meta)
nat_meta_s=nat_meta %>% dplyr::select ('ID','sample_barcode')
nat_meta_s=nat_meta_s[nat_meta_s$ID!='NA',]
nat_meta_s$ID=gsub('PT_','PT_C',nat_meta_s$ID)

colnames(nat_meta_s)[1]='cell_group'

final_table=rbind(nat_meta_s,tumor_map)
write.table(final_table,paste('TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/RNA_clusters',
'_Tumor_PT_clusters.snATAC.20210917.tsv',sep=''),sep='\t',row.names=F,quote=F)