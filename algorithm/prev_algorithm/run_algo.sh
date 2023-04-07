#!/bin/bash
#SBATCH --job-name=longitudinalAlgo
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=16GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./output0.out
#SBATCH --error=./error0.out
time python longitudinalTree.py -input ../Clustering/bnpc/gt_input_ForTree/default_rep5.Gprime.csv -cc ../Clustering/bnpc/sim_input/default/rep5/default_rep5.SNVcell.csv -celltp ../Clustering/bnpc/sim_input/default/rep5/default_rep5.SNVcell.csv -sim true
time python longitudinalTree.py -input ../Clustering/bnpc/sim_input_ForTree/default_rep1.Gprime.csv -cg ../Clustering/bnpc/sim_consensus_genotype/default_rep1.tsv -cc ../Clustering/bnpc/sim_input_ForTree/default_rep1.npy -D ../simulation/default/rep1/default_rep1.D.csv -celltp ../simulation/default/rep1/default_rep1.SNVcell.csv -alpha 0.01 -beta 0.1 -sim false

# Run the following for multiple runs
time python longitudinalTree.py -input ../Clustering/bnpc/sim_input_ForTree/default_mr_rep5.Gprime.csv -cg ../Clustering/bnpc/sim_output/default_mr/rep5/best/default_mr_rep5.G.tsv -cc ../Clustering/bnpc/sim_input_ForTree/default_mr_rep5.npy -D ../simulation/default/rep5/default_rep5.D.csv -celltp ../simulation/default/rep5/default_rep5.SNVcell.csv -alpha 0.01 -beta 0.1 -sim false -op_tree default_mr_op/rep5_tree.csv 
