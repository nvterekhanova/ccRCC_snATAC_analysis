system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")
	       
###For some samples (with many cells >6K) python-package used in chromVar doesn't work properly; Need to use this 2 commands:
###export OMP_NUM_THREADS=1
###export USE_SIMPLE_THREADED_LEVEL3=1
###export OPENBLAS_NUM_THREADS=1

library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


ATAC=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.',
'chromvar.cicero.v3.20210725.rds',sep=''))


###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type


#####################################
####Calculate score  per cell group##
#####################################
###Also separate NAT by Cluster
nat=readRDS(paste('../../10.Cell_of_origin/Map_snRNA_clusters.20210804/NAT_PT_LOH_cellsOnly.20210804.',
'v1.GeneActivity.rds.gz',sep=''))
nat_meta=nat@meta.data
cl=table(nat_meta$RNA_seurat_clusters[nat_meta$cell_type=='PT'])
pt_clusters_rna=names(cl[cl>200])
nat_meta$ID=paste(nat_meta$cell_type,nat_meta$RNA_seurat_clusters,sep='_')
nat_meta$ID=ifelse(nat_meta$RNA_seurat_clusters %in% pt_clusters_rna & nat_meta$cell_type=='PT',
nat_meta$ID, 'NA')


orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
nat_meta=nat_meta[orig_1$individual_barcode,]
ATAC$PT_cluster=nat_meta$ID

ATAC$test=ifelse(ATAC$cell_type!='PT',paste(ATAC$Piece_ID,ATAC$cell_type,sep='_'),ATAC$PT_cluster)

###Instead calculate for all cell types:
cell_types=unique(ATAC$test)


DefaultAssay(ATAC) <- 'chromvar'

chromv= GetAssayData(object = ATAC)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]


ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

for (cell_type in cell_types){
    if(length(rownames(ann_col1)[ann_col1$cell_types==cell_type])>=50){
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_type]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
	   mean_score=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           stat=cbind(cell_type,rownames(res)[motif],mean_score)
	   all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
 	all_wilcoxon_stat$mean_score=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score)))
	all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$mean_score),]
        final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
	print(cell_type)
}
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[2]='TF_Name'

write.table(final_wilcoxon_stat,paste("out/Motif_score_perCell_group.AllCellTypes.PTbyCluster.",
"20210804.tsv",sep=''),quote=FALSE,sep="\t",row.names=FALSE)
