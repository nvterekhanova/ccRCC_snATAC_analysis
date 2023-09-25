####2021-07-25:
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Signac)
library(Seurat)
library(reshape2)
library(ggfortify)
library(monocle3)
library(SeuratWrappers)


setwd('/Users/nadezhdaterekhanova/Desktop/Projects/ATAC/CCRCC/Analysis/Updated_pipeline/6.DA_motifs/For_cell_of_origin.20210721/Monocle_pseudotime.20210725/')

pt_nat=readRDS('../../../3.Merged_snATAC/Merge.v.20210707/NAT_PT_cellsOnly.20210721.v1.rds.gz')
PT=as.cell_data_set(pt_nat)
PT= cluster_cells(cds = PT, reduction_method = "UMAP")
PT<-learn_graph(PT,use_partition=T)
Idents(pt_nat)=pt_nat$seurat_clusters

#cluster_8=WhichCells(pt_nat, idents = 8)
#PT=order_cells(PT,reduction_method = "UMAP",root_cells=cluster_8)
#lets select cluster 8 manually instead
PT=order_cells(PT,reduction_method = "UMAP")
x=plot_cells(
cds = PT,color_cells_by = "pseudotime",
show_trajectory_graph = TRUE)


pdf('PT_cells_only_root_cluster8.Monocle3.20210725.pdf',width=6,height=5)
print(x)
dev.off()
