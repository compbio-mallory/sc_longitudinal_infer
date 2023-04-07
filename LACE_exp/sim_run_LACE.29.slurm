#!/bin/bash
#SBATCH --job-name=LACE29
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp4/rep5/output29.out
#SBATCH --error=sim_output/bp4/rep5/error29.out
module load R
Rscript run_LACE.R ./sim_input/bp4/rep5/bp4_rep5.RData sim_output/bp4/rep5/bp4_rep5_inf.rdata sim_output/bp4/rep5/bp4_rep5_lgTree.xml
