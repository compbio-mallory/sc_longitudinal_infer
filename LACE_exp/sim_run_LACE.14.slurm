#!/bin/bash
#SBATCH --job-name=LACE14
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap05/rep5/output14.out
#SBATCH --error=sim_output/ap05/rep5/error14.out
module load R
Rscript run_LACE.R ./sim_input/ap05/rep5/ap05_rep5.RData sim_output/ap05/rep5/ap05_rep5_inf.rdata sim_output/ap05/rep5/ap05_rep5_lgTree.xml
