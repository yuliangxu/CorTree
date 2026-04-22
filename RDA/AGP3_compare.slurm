#!/bin/bash
#SBATCH --job-name=AGP3cmp
#SBATCH --time=30:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40g
#SBATCH --cpus-per-task=1
#SBATCH --output=/hpc/home/yx306/CorTree/CorTree/out/out/AGP3_compare-%j.out
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab

module load R/4.4.0
R CMD BATCH --no-save --no-restore ./RDA/AGP3_compare.R ./out/AGP3_compare.Rout
