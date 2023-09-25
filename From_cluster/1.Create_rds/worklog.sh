export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

###provide samples:
for sample in Sample_1 Sample_2
do 
cellranger_output_directory=/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Cell_Ranger/outputs
macs2_path=/home/nvterekhanova/anaconda3/envs/r-environment/bin/macs2
Rscript atac_pipeline_v3.0.R -s $sample -a $cellranger_output_directory -m $macs2_path --prf_min 1000 --pct_min 15 --ns_max 5 --pc_first 2
done   
