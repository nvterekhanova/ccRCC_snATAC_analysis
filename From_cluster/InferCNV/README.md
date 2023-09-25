# Scripts for CNV analysis using **inferCNV**.

  * Reference tutorial: https://github.com/broadinstitute/inferCNV/wiki

  * Requires two input files:

   + Gene count matrix (raw counts).

   + Cell type annotation (Ref, Obs; or others).

## Steps:

1. ```Prepare_inputs/Prepare_inputs.R``` -- Prepare input gene count matrix.

  * Uses gene activity assay.

  * Need to use slot ```counts```, with not yet normalized data (was used differently before!) -- check tutorial above.

2. ```Script/run_infer_CNVs.sh``` -- General script for running inferCNV pipeline.


3. ```run_infer_CNVs_SampleID.sh``` -- Script for running inferCNV pipeline for a particular sample.