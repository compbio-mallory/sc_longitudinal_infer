#!/bin/bash
#SBATCH --job-name=LACE91
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u3/rep2/output91.out
#SBATCH --error=sim_output/u3/rep2/error91.out
module load R
Rscript run_LACE.R ./sim_input/u3/rep2/u3_rep2.RData sim_output/u3/rep2/u3_rep2_inf.rdata sim_output/u3/rep2/u3_rep2_lgTree.xml
