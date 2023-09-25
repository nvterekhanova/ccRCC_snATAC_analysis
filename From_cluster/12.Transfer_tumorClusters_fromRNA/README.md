# Tumor heterogeneity analysis.

  * Check also here ```PanCancer_snATAC/From_cluster/8.Tumor_heterogeneity/```.

## Analysis steps:

1. ```Subset_TumorOnly_reCluster.20210913.R``` -- Subset tumor cells for each sample and further re-clustering, dim. reduction.


2. ```TumorOnly_Integrate_withRNA.20210914.R``` -- Transfer snRNA clusters labels to tumor cells in snATAC objects.

  * It adds predictions to snATAC objects meta.data and saves them in ```TumorOnly_RNAClusters_mapped/metadata```. RDS are not saved, because annotation in the meta.data is enough.

  * Need to have snRNA clusters manual annotation table (done by Yige).


3. ```Re_assign_RNA_clusters.20210915.R``` -- Re-assigning of cluster ids using predicted annotation from label transfer (previous step).

  * It goes over all seurat clusters in ATAC, and checks if there are >50% of cells annotated with any cluster id. If yes, the cluster will be annotated with the respective cluster label.