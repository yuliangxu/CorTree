#!/bin/bash
#SBATCH --job-name=NRF1Sensi
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=20g
#SBATCH --cpus-per-task=1
#SBATCH --output=/hpc/home/yx306/CorTree/out/out/NRF1_sensi-%A_%a.out
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab
#SBATCH --array=1-3

module load R/4.4.0
Rscript ./RDA/NRF1_sensi.R "${SLURM_ARRAY_TASK_ID}" > ./out/NRF1_sensi${SLURM_ARRAY_TASK_ID}.Rout 2>&1
