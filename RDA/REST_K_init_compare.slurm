#!/bin/bash
#SBATCH --job-name=RESTKInit
#SBATCH --time=5:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=20g
#SBATCH --cpus-per-task=1
#SBATCH --output=/hpc/home/yx306/CorTree/out/out/REST_K_init_compare-%A_%a.out
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab
#SBATCH --array=1-4

module load R/4.4.0
Rscript ./RDA/REST_K_init_compare.R "${SLURM_ARRAY_TASK_ID}" > ./out/REST_K_init_compare${SLURM_ARRAY_TASK_ID}.Rout 2>&1
