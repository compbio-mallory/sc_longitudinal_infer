#!/bin/bash
#SBATCH --job-name=LACE28
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp4/rep4/output28.out
#SBATCH --error=sim_output/bp4/rep4/error28.out
module load R
Rscript run_LACE.R ./sim_input/bp4/rep4/bp4_rep4.RData sim_output/bp4/rep4/bp4_rep4_inf.rdata sim_output/bp4/rep4/bp4_rep4_lgTree.xml
