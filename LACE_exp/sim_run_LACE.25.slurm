#!/bin/bash
#SBATCH --job-name=LACE25
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp4/rep1/output25.out
#SBATCH --error=sim_output/bp4/rep1/error25.out
module load R
Rscript run_LACE.R ./sim_input/bp4/rep1/bp4_rep1.RData sim_output/bp4/rep1/bp4_rep1_inf.rdata sim_output/bp4/rep1/bp4_rep1_lgTree.xml
