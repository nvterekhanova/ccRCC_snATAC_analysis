# Scripts for performing TF motif enrichment analysis on bulk data.

  * Want to compare our results from SN study with the results from bulk ATAC study:

   + https://www.science.org/doi/10.1126/science.aav1898.

## Analyses:

### ```2021-04-09_ccRCC_peaks_ATAC``` -- This folder contains script for running DAP analysis on bulk.

#### Notes:    

  * ```2021-04-09_ccRCC_peaks_ATAC/Peaks_KIRC.20210409.R``` -- Calculate KIRC speicific DAPs.


  * For the above script need to download ```TCGA-ATAC_PanCan_Log2Norm_Counts.rds``` from here:

   + https://gdc.cancer.gov/about-data/publications/ATACseq-AWG.


  * The analysis is based on workshop here:

   + https://benbermanlab.com/assets/code/Workshop%20for%20ATAC-seq%20analysis.html

  * We calculate significantly UP peaks in comparisons below:

   1) KIRC vs all other cancers

   2) KIRC vs KIRP.

  * **Then take overlap of those two peaks lists**. 


  * Other steps were later updated, because ccRCC-specific peaks were updated, see below for latest scripts.


### ```2021-10-08_ccRCCpeaks_atac``` -- Version used right now, and for the Suppl. Fig. in the paper. Need to use selected DAPs from the above step.

#### Notes:

  * Need to use for it TFmotifView (online): http://bardet.u-strasbg.fr/tfmotifview/.

   + Use ref hg38 and paste all your peaks, that were identified as DAPs, into the window.

   + Perform two runs:

   1) without selection -- will run across all JASPAR motifs, but will use just one motif from each cluster of similar motifs, as stated in their *help message*. So, in the results of this run it won't have results for all TF motifs, that we found in SN analysis.

   2) we manually select several motifs of interest (the ones that were found to be ccRCC-specific when comparing with normal PT cells).

   + After results are ready, we download both tables, and can coombine them safely -- p-values are not affected.
    
    
  * Reference: https://academic.oup.com/nar/article/48/W1/W208/5824151.


3. ```11.Compare_with_Bulk``` -- First verson, not used right now.


