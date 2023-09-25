# Here are the scripts for merging snATAC samples for ccRCC project.

1. ```Merge_ATAC_samples_Recallpeaks.v.1.0.R``` -- merging samples using small number of peaks.

2. ```Overlap_peaks_acrossSamples.20210706.R``` -- making set of non-overlapping peaks, and quantifying coverage across those peaks for the merged object from the step 1.

3. ```RunChromVar.V2.R``` -- run chromVar to generate TF activity scores based on their motif accessibilites.