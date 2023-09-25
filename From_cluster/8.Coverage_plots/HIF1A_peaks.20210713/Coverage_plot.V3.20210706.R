#system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

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
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyverse)


markers=read_delim('ccRCC_markers.Surface.20210702.v1.csv',delim=',')
markers=as.data.frame(markers)
ATAC=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'20210707.rds',sep=''))
Idents(ATAC)=ATAC$Piece_ID

###Use the latest cell type annotation in "cell_type_manual_5" meta.data field:
cell_t=read.table('../../Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv',sep='\t',header=TRUE)

cell_t$individual_barcode=rownames(cell_t)

orig_1=as.data.frame(ATAC$dataset)
orig_1$individual_barcode=row.names(orig_1)
cell_t=cell_t[orig_1$individual_barcode,]
ATAC$cell_type=cell_t$cell_type


atac=ATAC
ATAC=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% c('')) | cell_type=='PT' & Piece_ID %in%
c(''))
 

meta=read.table('../Sample_categories.20210503.txt',sep='\t',header=T)
meta=meta[meta$Aliquot.WU %in% unique(ATAC$Piece_ID),]
bap1_s=meta$Aliquot.WU[meta$Category %in% c('BAP1-mutant','Both mutated')]
non_bap1_s=meta$Aliquot.WU[meta$Category %in% c('PBRM1-mutant','Non-mutant')]

####annotate all peaks:
peaks_1=StringToGRanges(rownames(ATAC@assays$peaksMACS2), sep = c("-", "-"))
 ###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks=anno
peaks$peak=paste(peaks$seqnames,peaks$start,peaks$end,sep='-')
peaks_sel=peaks[peaks$SYMBOL=='HIF1A' & peaks$annotation=='Promoter' & !is.na(peaks$SYMBOL),]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(ATAC) <- annotations

for (i in 1:nrow(peaks_sel)){
peak=peaks_sel$peak[i]
chr=gsub('(.*)-(.*)-(.*)','\\1',peak)
st=as.numeric(gsub('(.*)-(.*)-(.*)','\\2',peak))
en=as.numeric(gsub('(.*)-(.*)-(.*)','\\3',peak))
new_st=st-1000
new_en=en+1000
new_peak=paste(chr,new_st,new_en,sep='-')
gene=peaks_sel$SYMBOL[i]
p=CoveragePlot(
  object = ATAC,
  region = new_peak,
  annotation = TRUE,
  peaks = TRUE
)

pdf(paste("out/",gene,"_",new_peak,"Tumor_PT_cells.pdf",sep=""),width=8,height=8,useDingbats=FALSE)
print(p)
dev.off()
print(paste(new_peak,i,gene,sep=' '))
}


##############################
#####plot by Sample###########
##############################


unique_p=unique(peaks$Peak)
for (i in 1:length(unique_p)){
peak=unique_p[i]
chr=gsub('(.*)-(.*)-(.*)','\\1',peak)
st=as.numeric(gsub('(.*)-(.*)-(.*)','\\2',peak))
en=as.numeric(gsub('(.*)-(.*)-(.*)','\\3',peak))
new_st=st-1000
new_en=en+1000
new_peak=paste(chr,new_st,new_en,sep='-')

peaks_g=peaks[peaks$Peak==peak,]
tfs=c('ARNT::HIF1A','EGR1','HIF1A','HSF2','KLF14','KLF15','MXI1','NEUROD1','NEUROG2(var.2)',
'NFKB1','NFKB2','NRF1','RBPJ','SP1','SP3','SP9','SREBF2','TCFL5','ZBTB14','ZNF148','ZNF75D')
peaks_g=peaks_g[peaks_g$TF_name %in% tfs,]
gene=unique(peaks_g$Gene)
for (tf in peaks_g$TF_name){
peaks_tf=peaks_g[peaks_g$TF_name ==tf,]
if (nrow(peaks_tf)>0){
p=CoveragePlot(
  object = ATAC,
  region = new_peak,
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(peaks_tf$motif_coords, sep = c("-", "-")),
  ranges.title = paste(tf,"-motifs",sep=''),
  links=FALSE
)
dir.create(paste('VEGFA_peaks.BySample.20210520/',tf,sep=''))
pdf(paste("VEGFA_peaks.BySample.20210520/",tf,"/",gene,"_",tf,"_",new_peak,".pdf",
sep=""),width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()
}
print(paste(new_peak,i,gene,sep=' '))
}
}
#########


####Try to make CICERO PLOT:
####Links are absent?
conns=read.table('../7.CICERO/out/26_ccRCC_snATAC_CICERO.tsv',sep='\t',header=TRUE)
ccans=read.table('../7.CICERO/out/26_ccRCC_snATAC_CICERO.CCAN.tsv',sep='\t',header=TRUE)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(ATAC) <- links

####To many links... try different way:
#conns_1=conns[conns$Peak1 %in% unique_p | conns$Peak2 %in%]
conns_1=conns[conns$coaccess>0.5,]
conns_1=conns_1[!is.na(conns_1$Peak1),]
ccans_1=ccans[ccans$Peak %in% c(conns_1$Peak1,conns_1$Peak2),]
links_1<- ConnectionsToLinks(conns = conns_1, ccans = ccans_1)
Links(ATAC)<-links_1


p=CoveragePlot(
  object = ATAC,
  region ='chr6-43766858-44085257',
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(unique_p, sep = c("-", "-")),
  ranges.title = paste("VEGFA_peaks",sep=''),
  links=TRUE
)
pdf('VEGFA_peaks.BySample.20210520/Links_VEGFA.pdf',width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()
pdf('VEGFA_peaks.20210520/Links_Tumor_PT_VEGFA.pdf',width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()



p=CoveragePlot(
  object = ATAC,
  region ='chr3-52397000-52417000',
  annotation = TRUE,
  peaks = TRUE,
  links=FALSE
)
pdf('out/BAP1_coverage.pdf',width=8,height=9,useDingbats=FALSE)
print(p)
dev.off()




add=PeakPlot(
  ATAC,
  region=new_peak,
  color = "dimgrey",
  peaks=StringToGRanges(peaks_tf$motif_coordsc, sep = c("-", "-"))
)

x=CombineTracks(
  plotlist = list(p, add))


CoveragePlot(
  object = ATAC,
  region = "chr10-129465872-129468372",
#  features = "MGMT",
  annotation = TRUE,
  peaks = TRUE,
#  tile = TRUE,
  links = TRUE
)

#DotPlot(object=ATAC,features='MAPK8'

region='chr10-129466872-129467372'
peak_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey"
)
probe_all_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes
)
probe_bady_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes_bady
)
probe_ww_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes_ww
)
probe_more_plot=PeakPlot(
  ATAC,
  region=region,
  color = "dimgrey",
  peaks=probes_more
)
gene_plot <- AnnotationPlot(
  object = ATAC,
  region = region
)

cov_plot<-CoveragePlot(
  object = ATAC,
  region = region,
  annotation = FALSE,
  peaks = FALSE
)

CombineTracks(
  plotlist = list(cov_plot, peak_plot, probe_bady_plot,probe_ww_plot,probe_more_plot, gene_plot)
)


####V2:
CombineTracks(
  plotlist = list(cov_plot, peak_plot, probe_bady_plot,probe_ww_plot,probe_more_plot,probe_all_plot, gene_plot)
)
