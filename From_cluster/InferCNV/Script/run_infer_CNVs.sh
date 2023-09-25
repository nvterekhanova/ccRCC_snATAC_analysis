#Modified version using Yige's script

#!/bin/bash
#####Important: need to add to the docker run command this "--user $(id -u):$(id -g)"; otherwise doesn't working for me in the ccRCC_scratch (Nadya)

## specify the parent directory, under which input files and output directories are placed
### CHANGE THIS to /diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/ so the docker container can access everything under this directory, including reading inputs and make outputs
dir_run=/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/InferCNV/

## specify the docker image, this is an environment with all the dependent R packages installed
docker_image_name=singlecellportal/infercnv

## specify the inferCNV running mode
analysis_mode=subclusters

## specify the path to the raw count matrix
### specify the aliquot id to enable getting the path to the raw count matrix
### CHANGE THIS for each aliquot
snRNA_aliquot_id=CPT0025880013
### specify the directory with all the raw count matrix
dir_raw_count_file=/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/InferCNV/inputs/raw_counts_matrix/
### specify the path to the raw count matrix for this aliquot
path_raw_count_file=${dir_raw_count_file}${snRNA_aliquot_id}.RNA_Count.tsv

## specify the path to the annotation file
### specify the run id for this run because the annotation file might change in the future
run_id=20200405.v2
### specify the directory with all the annotation files
dir_annotation_file=/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/InferCNV/inputs/annotations_file/${run_id}/
### specify the path with the annotation file for this aliquot
path_annotation_file=${dir_annotation_file}${snRNA_aliquot_id}.Barcode_Annotation.txt

## specify the path to the gene order file, this file contains the position information of each gene
path_gene_order_file=/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt

## specify the cutoff for the average gene “expression” level to filter genes
cutoff=0.04

## specify the directory to the outputs for this aliquot
### specify the directory to the outputs for all the runs
dir_infercnv_outputs=/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/InferCNV/outputs/
### specify the directory to the outputs for this run
dir_infercnv_outputs_by_run=${dir_infercnv_outputs}${run_id}/
mkdir -p ${dir_infercnv_outputs_by_run}
### specify the directory to the outputs for this run and for this aliquot
dir_output=${dir_infercnv_outputs_by_run}${snRNA_aliquot_id}/
mkdir -p ${dir_output}

## specify the label for the reference cells in the annotation file
### for reference cells with different cell types, please combine them here by “,”, for example “Endothelial cells,Macrophages,CD4+ T-cells"
ref_group_names=DC,Fibroblasts,Macrophages_M1,CD8+_T-cells,Macrophages_M2,Endothelial_cells

## specify the path to the log file, each run has a unique timestamp
path_log_file=${dir_output}${snRNA_aliquot_id}.$(date +%Y%m%d%H%M%S).log

## below is the main running code
### the following line specify 
docker run --user $(id -u):$(id -g) -v ${dir_run}:${dir_run} ${docker_image_name} inferCNV.R \
        --analysis_mode=${analysis_mode} \
        --raw_counts_matrix=${path_raw_count_file} \
        --annotations_file=${path_annotation_file} \
        --gene_order_file=${path_gene_order_file} \
        --cutoff=${cutoff} \
        --out_dir=${dir_output} \
        --cluster_by_groups \
        --denoise \
        --HMM \
        --num_threads=20 \
        --ref_group_names=${ref_group_names} &> ${path_log_file}&
