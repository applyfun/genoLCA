#!/bin/sh
#SBATCH --time=01:00:00
#SBATCH -p shared,brc
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --ntasks=8

module load apps/R/3.6.0

cd ~/brc_scratch/scripts/genotype_lca/

R CMD BATCH --no-environ --no-save "simulate_poLCA_data.r" "outfile_simulation_lca_data.txt"
