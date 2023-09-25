#script for the marker-plotting: "../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/Check_PT_segmental.markers.R"

library(Signac)
library(Seurat)

ATAC=readRDS(paste('../../3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.',
'cicero.v3.20210725.rds',sep=''))

peak.data <- GetAssayData(object = ATAC, assay = 'peaksMACS2', slot = "counts")
total_fragments_cell <- ATAC$passed_filters
peak.counts <- colSums(x = peak.data)
frip <- peak.counts / total_fragments_cell
ATAC <- AddMetaData(object = ATAC, metadata = frip, col.name = 'frip_500MACS2')
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_500MACS2')

barcode_map=read.table(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/',
'12.Transfer_tumorClusters_fromRNA/TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/',
'RNA_clusters_Tumor_PT_clusters.snATAC.20210917.tsv',sep=''),sep='\t',header=T)
rownames(barcode_map)=barcode_map$sample_barcode

atac=ATAC
ATAC=subset(atac, cells=barcode_map$sample_barcode)

orig_1=as.data.frame(ATAC$orig.ident)
barcode_map=barcode_map[rownames(orig_1),]
ATAC$cell_group_ID=barcode_map$cell_group
Idents(ATAC)=ATAC$cell_group_ID

DefaultAssay(ATAC)<-'ATACGeneActivity'

###Calculate average GeneActivity for selected markers EMT-markers:
markers=read.table('data/EMT_vs_selectedEpithelialClusters.DEG_overlap_DAP.Consistent.20210928.v1.tsv',header=T,sep='\t')

aliquot.averages <- AverageExpression(ATAC, assays = 'ATACGeneActivity', slot ='data',features=markers$Gene)
#Warning: The following 38 features were not found in the ATACGeneActivity assay: ARHGAP24, ATP5MC2, AUTS2, CCDC171, DRAIC, EGLN3, EIPR1, EMX2OS, FOXP2, FRG1-DT, GLIS3, GPHN, GRAMD2B, GUCY1B1, HDHD5, JPT2, LINC00909, LINC01320, LINC01508, LINC02532, LRBA, LUCAT1, MAST4, MIATNB, MMP24OS, MYOCOS, NEAT1, PDE10A, RETREG1, SAP30L-AS1, SMC5-AS1, SNHG15, TENT4A, TRMT9B, TUT7, VPS35L, VXN, XIST

file2write <- paste0("out/AverageGeneActivity.EMT_vs_selectedEpithelialClusters.DEG_overlap_DAP.Consistent.20210928.tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)

###Calculate average GeneActivity for selected markers ccRCC cell surface-markers:
markers=read_delim('data/ccRCC_markers.Surface.20210824.v1.tsv',delim='\t')
markers=as.data.frame(markers)

aliquot.averages <- AverageExpression(ATAC, assays = 'ATACGeneActivity', slot ='data',features=markers$Gene)
#Warning: The following 3 features were not found in the ATACGeneActivity assay: DPP6, EPHA6, PLCB1

file2write <- paste0("out/AverageGeneActivity.ccRCC_markers.Surface.20210928.tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)


###Also make a dotplot, just to check:
#x1=DotPlot(ATAC,features=c('SLC5A12', 'SLC13A3', 'SLC5A2', 'SLC3A1', 'SLC38A3', 'SLC16A9'))+ RotatedAxis()

#pdf("out/S1_S2_S3_marker_plots.v2.20210917.pdf",width=17,height=10)
#print(x1)
#dev.off()
