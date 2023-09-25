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
library(SeuratWrappers)
library(patchwork)
library(cicero)
library(EnsDb.Hsapiens.v86)

###we also ran the GeneActivity on the object
ATAC=readRDS(paste('../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'20210707.rds',sep=''))

DefaultAssay(ATAC)<-'peaksMACS2'
ATAC.cds <- as.cell_data_set(x = ATAC)
ATAC.cicero <- make_cicero_cds(ATAC.cds, reduced_coordinates = reducedDims(ATAC.cds)$UMAP)

#Overlap QC metrics:
#Cells per bin: 50
#Maximum shared cells bin-bin: 44
#Mean shared cells bin-bin: 0.0112756486209922
#Median shared cells bin-bin: 0

# get the chromosome sizes from the Seurat object
genome = seqlengths(BSgenome.Hsapiens.UCSC.hg38)

# use chromosome 1 to save some time -- we don't do it, because need for whole genome
# omit this step to run on the whole genome
#genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)
chrs=paste('chr',c(1:22,'X','Y'),sep='')
genome.df=genome.df[genome.df$chr %in% chrs,]

# run cicero
conns <- run_cicero(ATAC.cicero, genomic_coords = genome.df, sample_num = 100)

write.table(conns, "out/28_ccRCC_snATAC_CICERO.tsv",sep='\t',quote=FALSE,row.names=FALSE)

ccans <- generate_ccans(conns)
#"Coaccessibility cutoff used: 0.3"

write.table(ccans, "out/28_ccRCC_snATAC_CICERO.CCAN.tsv",sep='\t',quote=FALSE,row.names=FALSE)

conns_1=conns[conns$coaccess>0.25,]
conns_1=conns_1[!is.na(conns_1$coaccess),]

write.table(conns_1, "out/28_ccRCC_snATAC_CICERO.0.25_cutoff.tsv",sep='\t',quote=FALSE,row.names=FALSE)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(ATAC) <- links

saveRDS(ATAC,paste("../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.",
"chromvar.cicero.v3.20210725.rds",sep=''),compress=T)



####Now annotate peaks:
library(Signac)
library(Seurat)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(reshape)
library(motifmatchr)

tab=read.table('out/28_ccRCC_snATAC_CICERO.0.25_cutoff.tsv',sep='\t',header=TRUE)
peaks_1=StringToGRanges(tab$Peak1, sep = c("-", "-"))
###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")

anno=as.data.frame(peakAnno)
peaks=tab
peaks$Gene=anno$SYMBOL
peaks$Type=anno$annotation
peaks$geneId=anno$geneId
peaks$motif_distanceToTSS=anno$distanceToTSS
peaks$geneStart=anno$geneStart
peaks$geneEnd=anno$geneEnd
write.table(peaks, 'out/28_ccRCC_snATAC_CICERO.0.25_cutoff.Annotated_peaks.tsv',
sep='\t',quote=FALSE,row.names=FALSE)


#######################
###ENDS HERE FOR NOW###
#######################









