#!/bin/bash
#SBATCH --job-name=tree_sim1_cortree
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=10g
#SBATCH --cpus-per-task=1
#SBATCH --output=/hpc/home/yx306/CorTree/CorTree/out/out/sim1-compare-%j.out
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab
#SBATCH --array=1-300

# srun -c 4 --pty bash -i
module load R/4.4.0
R CMD BATCH --no-save --no-restore ./test/Sim1_clus_cortree.R  ./out/out/Sim1_clus${SLURM_ARRAY_TASK_ID}.Rout


