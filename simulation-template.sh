#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH -p shared,brc
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --ntasks=8

module load apps/R/3.6.0

cd ~/brc_scratch/scripts/genotype_lca/

R CMD BATCH --no-environ --no-save '--args IDX' run_LCA.r outfile_lca_modelling_IDX.txt
