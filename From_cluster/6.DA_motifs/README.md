# DAM analysis.

  * Calculate DAMs between CCRCC and PT

   + ```Calculate_Motif_score.R``` is used to calculate DAMs, and also average TF motif score per cell group. Also see example here: ```HTAN_BRCA/From_cluster/5.DA_motifs/Subtype_DAMs.20211203/```. Though, it was calculated in a manner of each tumor sample vs pooled PT cell together, and then consisently UP/DOWN TF motifs were selected.


  * Calculate average motif score per celll group.

```Scores_whenSeparatingPT_byCluster_snRNA_Clusters/Calculate_Motif_score.PT_byCluster.R``` -- this script calculates average score when separating PT by cluster.