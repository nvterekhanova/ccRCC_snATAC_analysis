###We want to re-size peaaks relative to summits +/-250bp (the same as in the ChromVar and CGA ATAC paper)

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 100 * 1024 ^ 2)

####Subsetting of the PT and LOH-clusters


######################################################
#Subset furtther PT/LOH-cells, and try to re-cluster##
######################################################

####UPD:Neeed to subset large object instead by cell ids, because can't calc GeneActivity on NAT-only object
#Loop of Henle            PT
#         6248         12335
orig_1=as.data.frame(nat$Piece_ID)
orig_1$Barcode=rownames(orig_1)
colnames(orig_1)[1]='Piece_ID'
write.table(orig_1,"LOH_PT_cells_barcode_ids.tsv",sep='\t',row.names=T,quote=F)

###Subset NATs:
atac=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/',
'28_ccRCC_snATAC.selectedPeaks.chromvar.cicero.v3.20210725.rds',sep=''))

nat=subset(atac,cells=orig_1$Barcode)

###When you subset the cells the inverse document frequency value (IDF) used in TF-IDF will change, so I would recommend re-running RunTFIDF() and RunSVD() on the subset data. (https://github.com/timoast/signac/discussions/565)

nat <- FindTopFeatures(nat, min.cutoff = 'q0')
nat <- RunTFIDF(nat)
nat <- RunSVD(nat)

nat <- RunUMAP(object = nat, reduction = 'lsi', dims = 2:50)
nat <- FindNeighbors(object = nat, reduction = 'lsi', dims = 2:50)


###this helps:
options(future.globals.maxSize= 891289600)
nat <- FindClusters(object = nat, verbose = FALSE, algorithm = 3)

nat$Piece_ID_1=ifelse(nat$Piece_ID %in% c(),nat$Piece_ID,'Other') 
pdf("PT_LOH_cells_reClustered_Cluster.50PCs.pdf",height=7,width=16,useDingbats = F)
p1=DimPlot(nat,group.by='Piece_ID_1')
p2=DimPlot(nat,group.by='seurat_clusters',label=T)
p1+p2
dev.off() 
pdf("PT_LOH_cells_reClustered_Cell_type.50PCs.pdf",height=7,width=16,useDingbats = F)
p1=DimPlot(nat,group.by='Piece_ID_1')
p2=DimPlot(nat,group.by='cell_type',label=T)
p1+p2
dev.off() 

####Now try to integraate Ruiyang's object:
rna_nat=readRDS(paste('/diskmnt/Projects/Users/rliu/Projects/ccRCC/PT_LOH_subcluster_analysis/snRNA/',
'merged_ccRCC_34_samples/analysis/subcluster_analysis/PT/results/PT_LOH.rds',sep=''))


pdf("PT_LOH_cells_snRNA.pdf",height=7,width=16,useDingbats = F)
p1=DimPlot(rna_nat,group.by='Sample')
p2=DimPlot(rna_nat,group.by='cell_type',label=T)
p1+p2
dev.off() 
pdf("PT_LOH_cells_ByCluster.snRNA.pdf",height=7,width=16,useDingbats = F)
p1=DimPlot(rna_nat,group.by='Sample')
p2=DimPlot(rna_nat,group.by='seurat_clusters',label=T)
p1+p2
dev.off() 


####Try to integrate and do label transfer for two datasets:
DefaultAssay(nat) <- 'ATACGeneActivity'
ATAC=nat
RNA=rna_nat

ATAC<- NormalizeData(
  object = ATAC,
  assay = 'ATACGeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC$nCount_ATACGeneActivity)
)

DefaultAssay(RNA)='RNA'
RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA<-FindVariableFeatures(RNA)
RNA<-ScaleData(RNA)

transfer.anchors <- FindTransferAnchors(
  reference = RNA,
  query = ATAC,
  reduction = 'cca'
)

RNA$RNA_seurat_clusters=RNA$seurat_clusters
RNA$RNA_cell_type=RNA$cell_type

annotations=c('RNA_seurat_clusters','RNA_cell_type')

for (label in annotations){
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = factor(RNA@meta.data[,colnames(RNA@meta.data)== label]),
  weight.reduction = ATAC[['lsi']],
  dims=2:50
)

ATAC <- AddMetaData(object = ATAC, metadata = predicted.labels)
colnames(ATAC@meta.data)=gsub('predicted.id',label,colnames(ATAC@meta.data))
}


####Make Dimplots:
pdf("plots/RNA_clusters_PT_LOH_cells_reClustered_Cluster.50PCs.pdf",height=7,width=16,useDingbats = F)
p1=DimPlot(ATAC,group.by='Piece_ID_1')
p2=DimPlot(ATAC,group.by='RNA_seurat_clusters',label=T)
p1+p2
dev.off() 
pdf("plots/RNA_cell_type_PT_LOH_cells_reClustered_Cell_type.50PCs.pdf",height=7,width=16,useDingbats = F)
p1=DimPlot(ATAC,group.by='Piece_ID_1')
p2=DimPlot(ATAC,group.by='RNA_cell_type',label=T)
p1+p2
dev.off() 



saveRDS(ATAC,"NAT_PT_LOH_cellsOnly.20210804.v1.GeneActivity.rds.gz",compress=TRUE)

#######################
####ENDS HERE##########
#######################

