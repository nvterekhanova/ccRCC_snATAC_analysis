## Reference: https://rdrr.io/cran/Seurat/src/R/differential_expression.R#sym-LRDETest
####Importnt: load Seurat after	Signac

library(Signac)
library(Seurat)
library(lmtest)
library(future.apply)


plan("multiprocess", workers =5)
options(future.globals.maxSize = 15 * 1024^3)

###This is an RDS file with columns corresponding to peaks, and rows to barcodes
barcode=readRDS('data/Barcode2BAP1_vs_Nonmutants_PrefilteredPeak.CNV.20210712.v1.RDS')

###This is an RDS file with columns corresponding to peaks, and rows to cases 
#case=readRDS('data/Case2Peak.CNV.20210604.v1.RDS')


ATAC$test='Other'
ATAC$test=ifelse(ATAC$cell_type=="Tumor" & ATAC$Piece_ID %in% bap1_s,"BAP1_mutants",ATAC$test)
ATAC$test=ifelse(ATAC$cell_type=="Tumor" & ATAC$Piece_ID %in% non_bap1_s,"NOT_BAP1_mutants",ATAC$test)

###Now check how FindMarkers modules will work:
### name it as "object"
object=ATAC
## set ident.use.1
ident.use.1='BAP1_mutants'
## set ident.use.2
ident.use.2='NOT_BAP1_mutants'
## make sure set Idents for the object
Idents(object)=object$test

##assay
assay='peaksMACS2'

## input the CNV values per peak per cell
### column is named after the peak by rownames(object)
### row is named after the cells by colnames(object)
cnv_per_feature_df=barcode


# set pre-filtering parameters --------------------------------------------
## set min.pct
min.pct=0.1
## set min.diff.pct (lower to 0, since it filtered a lot of peaks)
min.diff.pct=0
## set logfc.threshold
logfc.threshold=0

# feature selection -------------------------------------------------------
fc.results <- FoldChange(
  object = object,
  slot = "data",
  ident.1 = ident.use.1,
  ident.2 = ident.use.2,
  base = 2
)

# feature selection (based on percentages)
alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
names(x = alpha.min) <- rownames(x = fc.results)
features <- names(x = which(x = alpha.min >= min.pct))
if (length(x = features) == 0) {
  warning("No features pass min.pct threshold; returning empty data.frame")
  return(fc.results[features, ])
}
alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
features <- names(
  x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
)
if (length(x = features) == 0) {
  warning("No features pass min.diff.pct threshold; returning empty data.frame")
  return(fc.results[features, ])
}
# feature selection (based on logFC)
#slot
slot='data'
only.pos=FALSE

if (slot != "scale.data") {
  total.diff <- fc.results[, 1] #first column is logFC
  names(total.diff) <- rownames(fc.results)
  features.diff <- if (only.pos) {
    names(x = which(x = total.diff >= logfc.threshold))
  } else {
    names(x = which(x = abs(x = total.diff) >= logfc.threshold))
  }
  features <- intersect(x = features, y = features.diff)
  if (length(x = features) == 0) {
    warning("No features pass logfc.threshold threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
}

filtered.peaks=fc.results[rownames(fc.results) %in% features,]
filtered.peaks$peak=rownames(filtered.peaks)
write.table(filtered.peaks,'out/BAP1_vsNon_mutants.Filtered_peaks_byMinPct.0.1.20210712.tsv',sep='\t',
row.names=F,quote=F)

## set "features" -- we do it later
features=colnames(barcode)


# set inputs for the below process --------------------------------------------------------------
cells.1 <- WhichCells(object = object, idents = ident.use.1)
cells.2 <- WhichCells(object = object, idents = ident.use.2)
data.use=object[[assay]]
data.use <- data.use[features, c(cells.1, cells.2)]

# process inputs further --------------------------------------------------
## prepare latent.vars data frame
latent.vars <- FetchData(
  object = object,
  vars = "peak_RF_500MACS2",
  cells = c(cells.1, cells.2)
)
## prepare group.info data frame
group.info <- data.frame(row.names = c(cells.1, cells.2))
group.info[cells.1, "group"] <- "Group1"
group.info[cells.2, "group"] <- "Group2"
group.info[, "group"] <- factor(x = group.info[, "group"])

## prepare data.use object
data.use <- data.use[, rownames(group.info), drop = FALSE]
latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
latent.vars <- cbind(latent.vars, cnv_per_feature_df[rownames(group.info),])
colnames(latent.vars)=gsub('-','_',colnames(latent.vars))
rownames(data.use)=gsub('-','_',rownames(data.use))

# run test ----------------------------------------------------------------
my.sapply <- ifelse(
  nbrOfWorkers() == 1,
#  test = verbose && nbrOfWorkers() == 1,
  yes = pbsapply,
  no = future_sapply
)
p_val <- my.sapply(
  X = rownames(x = data.use),
  FUN = function(x) {
    model.data <- cbind(GENE = data.use[x,], group.info, latent.vars)
    fmla <- as.formula(object = paste(
      "group ~ GENE +",
      paste(c(x, 'peak_RF_500MACS2'), collapse = "+")
    ))
    fmla2 <- as.formula(object = paste(
      "group ~",
      paste(c(x, 'peak_RF_500MACS2'), collapse = "+")
    ))
    model1 <- glm(formula = fmla, data = model.data, family = "binomial")
    model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
    lrtest <- lrtest(model1, model2)
    return(lrtest$Pr[2])
  }
)
to.return <- data.frame(p_val, row.names = rownames(data.use))
to.return$p_val=as.numeric(as.character(unlist(to.return$p_val)))
to.return$FDR=p.adjust(to.return$p_val,method='fdr')
to.return$p_adjust_bonf=p.adjust(to.return$p_val,method='bonferroni')

###Now merge with fc.results
to.return$peak=row.names(to.return)
fc.results$peak=row.names(fc.results)
to.return$peak=gsub('_','-',to.return$peak)
to.return=merge(to.return,fc.results,all.x=TRUE)


to.return$chr_peak=gsub('(.*)-.*-.*','\\1',to.return$peak)

write.table(to.return,"out/DA_peaks_BAP1mutants_vs_nonMutants_correctedbyCNV.20210713.tsv",
sep='\t',quote=FALSE,row.names=F)