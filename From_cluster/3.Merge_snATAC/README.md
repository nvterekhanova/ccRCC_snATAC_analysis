# Scripts used to create merged objects using all CCRCC cohort samples, including snATAC and multiome.

## Note:

   * The scripts in here were slightly updated since the time they were used for CCRCC project, for up-to-date scripts use (both are same):

   + ```HTAN_BRCA/From_cluster/3.Merge_snATAC/```

   + ```SC_analysis_scripts/3.Merge_objects/```


## Steps include:

1. ```Merge_ATAC_samples_Recallpeaks.v.1.0.R``` -- Make prelim cohort level objects using just small sets of peaks.

2. ```Overlap_peaks_acrossSamples.20210706.R```

   * Combine peaks from all samples, and then remove overlapping ones.

   * Quantify new set of peaks and make final merged object.

3. ```RunChromVar.V2.R``` -- Calculate motif-scores with chromVAR, using Seurat wrapper for it, script ```RunChromVar.V2.R```. Save the object with the added assays -- that will be the final object.
