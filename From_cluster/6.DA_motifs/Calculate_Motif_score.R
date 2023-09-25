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


ATAC=readRDS(paste('../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'20210707.rds',sep=''))


###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type

atac=ATAC
ATAC=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% c('C3L-00004-T1', 'C3L-00010-T1',
 'C3L-00079-T1', 'C3L-00088-T2', 'C3L-00088-T1',  'C3L-00448-T1',
'C3L-00583-T1', 'C3L-00610-T1', 'C3L-00790-T1', 'C3L-00908-T1', 'C3L-00917-T1',
'C3L-01287-T1',  'C3L-01313-T1', 'C3N-00317-T1','C3N-00733-T1',
'C3N-01200-T1','C3N-01213-T1','C3L-00416-T2','C3L-01302-T1','C3L-00096-T1', 'C3L-00026-T1',
'C3N-00242-T1', 'C3N-00437-T1', 'C3N-00495-T1')) | cell_type=='PT' & Piece_ID %in%
c('C3L-00088-N','C3N-01200-N','C3L-00079-N','C3N-00242-N'))


###Tumor-Normal comparison:

ATAC$test=ifelse(as.character(ATAC$Piece_ID) %in% c("C3N-01200-N","C3L-00088-N","C3L-00079-N","C3N-00242-N"),"PT",
paste(as.character(ATAC$Piece_ID),"Tumor",sep="_"))

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

tumor_piece_ids=unique(ATAC$Piece_ID)
tumor_piece_ids=tumor_piece_ids[!(tumor_piece_ids %in% c("C3L-00088-N","C3N-01200-N","C3L-00079-N","C3N-00242-N"))]
cell_types=paste(tumor_piece_ids,"Tumor",sep='_')


final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

for (cell_t1 in cell_types){
    print (cell_t1)
    cell_t2='PT'
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t2]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")


final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}


final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)	
colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')
#final_wilcoxon_stat$diff=final_wilcoxon_stat$mean_score2-final_wilcoxon_stat$mean_score1
#final_wilcoxon_stat=final_wilcoxon_stat[order(-final_wilcoxon_stat$diff),]

write.table(final_wilcoxon_stat,paste("out/Score_difference.Tumor_Normal_comparison.20210715.tsv",
sep=""),quote=FALSE,sep="\t",row.names=FALSE)


#####################################
####Calculate score  per cell group##
#####################################

ATAC$test=paste(ATAC$Piece_ID,ATAC$cell_type,sep='_')

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

write.table(final_wilcoxon_stat,paste("out/Motif_score_perCell_group.AllCellTypes.20210720.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)
