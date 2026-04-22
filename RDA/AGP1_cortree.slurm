#!/bin/bash
#SBATCH --job-name=AGP1
#SBATCH --time=30:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30g
#SBATCH --cpus-per-task=1
#SBATCH --output=/hpc/home/yx306/CorTree/CorTree/out/out/AGP1-%A_%a.out
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab
#SBATCH --array=1

# srun -c 4 --pty bash -i
module load R/4.4.0
R CMD BATCH --no-save --no-restore ./RDA/AGP1_cortree.R ./out/AGP1_cortree${SLURM_ARRAY_TASK_ID}.Rout