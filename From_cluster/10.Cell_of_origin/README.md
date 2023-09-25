# Analyses related to cell of origin.

## Notes:

1. We subtract normal PT and LOH cells from the merged object, and do re-clustering.

   * Merged object should contain already Gene Activities, otherwise would be return error when you try to calculate them on the subsetted object, issue that wasn't possible to fix quickly.


2. Next we use procedure of lable transfer, similar to the one used for tumor heterogeneity analysis, in order to annotate clusters in snATAC with the respective snRNA cluster labels. This is needed to have corresponding clusters between snATAC/snRNA.

  * Check here as well: ```PanCancer_snATAC/From_cluster/8.Tumor_heterogeneity/```.