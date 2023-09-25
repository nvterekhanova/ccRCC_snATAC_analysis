library(Signac)
library(Seurat)

tab=read.table('DEG_associated_Peaks.onlyPromoters_Enhancers.20210514.v1.tsv',sep='\t',header=TRUE)
enriched.motifs <- FindMotifs(
  object = ATAC,
  features = tab$Peak
)

enriched.motifs$motif_name=row.names(enriched.motifs)
write.table(enriched.motifs,paste("Enriched_motifs_inDEG_associated_Peaks.onlyPromoters_Enhancers.",
"20210514.v1.tsv",sep=''),sep='\t',quote=FALSE,row.names=FALSE)

write.table(tab,paste("DEG_associated_Peaks.onlyPromoters_Enhancers.",
"20210514.v1.tsv",sep=''),sep='\t',quote=FALSE,row.names=FALSE)