###For this analysis we don't need pt, only subset Tumor cells (same way as was done for snRNA)
#for annotation use code from /diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/6.DA_motifs/Scores_whenSeparatingPT_byCluster_snRNA_Clusters
###conda enironment r_4:
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 100 * 1024 ^ 2)



cell_t=read.table('../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

samples=c('')

for (sample in samples){
atac=readRDS(paste('../1.Create_rds/out/',sample,'/',sample,'_processed_atac.rds',sep=''))
cell_t_sample=cell_t[cell_t$dataset==sample,]

rownames(cell_t_sample)=gsub('(.*)_(.*)','\\2',rownames(cell_t_sample))
orig_1=as.data.frame(atac$orig.ident)
cell_t_sample=cell_t_sample[rownames(orig_1),]

atac$cell_type=cell_t_sample$cell_type

atac_s=subset(atac,cell_type %in% c('Tumor','EMT tumor cells'))

###Now re-cluster:
atac_s <- FindTopFeatures(atac_s, min.cutoff = 'q0')
atac_s <- RunTFIDF(atac_s)
atac_s <- RunSVD(atac_s)

atac_s <- RunUMAP(object = atac_s, reduction = 'lsi', dims = 2:30)
atac_s <- FindNeighbors(object = atac_s, reduction = 'lsi', dims = 2:30)


###this helps:
options(future.globals.maxSize= 891289600)
atac_s <- FindClusters(object = atac_s, verbose = FALSE, algorithm = 3)

plot1=DimPlot(atac_s, group.by='seurat_clusters',label=T)
plot2=DimPlot(atac_s, group.by='cell_type',label=T)
p=plot1+plot2

pdf(paste('TumorOnly.RDS.dimplots/',sample,'_TumorOnly.20210913.pdf',sep=''),width=12,height=5)
print(p)
dev.off()

saveRDS(atac_s,paste('TumorOnly.RDS/',sample,'_TumorOnly.20210913.rds',sep=''))
print(sample)
}
