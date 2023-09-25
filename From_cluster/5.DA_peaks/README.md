# Analyses scripts related to peak accesibility.

## List of analyses:

1. ```ccRCC_specific_peaks.20210711``` -- CCRCC specific peaks, DAP analysis between CCRCC tumor and normal prox. tubule (PT) cells.


2. ```Cell_typeSpecific_peaks.20210811``` -- Cell type specific DAPs. Comparison one cell type vs all others.


3. ```BAP1_da_peaks.20210608``` -- BAP1 specific DAPs. Comparison of BAP1-mutant samples vs non-BAP1/PBRM1 mutants, tumor cells only.

  * Two samples with both BAP1 and PBRM1 mutations considered as BAP1-mutants, because of larger pronounced effect of BAP1 (based on literature and observation).


4. ```PBRM1_da_peaks.20210615``` -- PBRM1 specific DAPs. Comparison of PBRM1-mutant samples vs non-BAP1/PBRM1 mutants, tumor cells only.

## Notes for analyses from 1, 3, 4:

  * We perform two runs of ```FindMarkers```: with CNV-normalization and without it, see details below.

## Steps and details:

  * Use for example ```BAP1_da_peaks.20210608/BAP1_specific_peaks_against_NonMutants/```.

1. ```Calculate_DA_peaks.20210622.R``` -- Run DAP analysis for each tumor sample from BAP1-mutants group vs pooled tumor cells from other group (non-BAP1/PBRM1 mutants).

2. ```FindMarkers_withCNV_correction.20210622.R``` -- Run DAP analysis for pooled tumor cells from BAP1-mutants group vs pooled tumor cells from other group (non-BAP1/PBRM1 mutants).

  * For this modified version of ```FindMarkers``` function is used. It includes additional covariate of CNV matrix (CNV per cell per peak)

   + CNV-matrix was calculated based on bulk data (by Yige), gene-level CNVs were used, and annotation of peaks based on nearest gene, e.g. promoter, intron, etc.

  * Further steps are performed on local PC, check ```ccRCC_SN_project/From_local/5.DA_peaks```.

   + The idea is tha we check for concordantly UP/DOWN peaks across results of comparisons from **step 1**, and then take only DAPs that were found significant in **step 2** as well (and the ones that are absent, if they were not present in the input set of peaks -- becuase CNVs were calculated not for all peaks (?) ). 


3. ```Peaks_coverage.CountsNorm.20211011.R``` -- The most updated version of average accessibility calculation.

  * Uses slot ```counts``` normalization instead of data, though results seem similar when use ```data``` slot.

  * ```Peaks_coverage.R``` -- Previous version of script using ```data``` slot.