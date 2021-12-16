# sc_longitudinal_infer

For clustering we used SCG. Here are the steps to run SCG:

1. conda activate scg
2. scg run_singlet_model --config_file examples/config.yaml --out_dir ../no_doublet

The files generated are cluster_posteriors.tsv and genotype_posteriors.tsv. Unzip the * tsv.gz to get the recent results. 
The cluster numbers can be changed in the config.yaml file. 

Simulated time points for each cell is obtained from SCG clusters (no. of clusters 5) by running this file: timepoints_from_cluster.py. 

Cells belonging to each clusters is obtained by running the file: cellToCluster.py. 

A json file is returned which contains the information and this file is used by 
equation2Value.py to get the information. (can change this later)

--------------------------------------------------------------------------------------

Simulated input for the longitudinal tree can be obtained by running the following command:

python input_for_tree.py -input ./examples/example_input_D.tsv -gp ./examples/example_genotype_posterior.tsv -cellCluster ./examples/example_cellToCluster.json -timepoints ./examples/example_cell_timepoints.json

Inputs required: an input D matrix, genotype_posterior file similar to output of SCG, cells assigned to clusters information in example_cellToCluster.json and the cell timepoints mentioned in example_cell_timepoints.json.

Output: a json file containing subclones arranged in timepoints. This output can be used as an input for the longitudinal tree algorithm.

---------------------------------------------------------------------------------------

The Equation 2 is implemented in the file equation2Value.py and it accepts the following parameters as an input: the input data matrix, two output files from SCG (genotype_posterior.tsv and cluster_posterior.tsv), the lambda coefficient value and number of clusters. This file can be executed with the following command:

python equation2Value.py -input ../scg/examples/snv.tsv -gp ../no_doublet/genotype_posteriors.tsv -cp ../no_doublet/cluster_posteriors.tsv -lambda_coeff 0.2 -subclones 40

Please replace it with corresponding path, file names and parameters.

---------------------------------------------------------------------------------

To run the Longitudinal Tree program use the following command: python longitudinalTree.py -input example_InputGprime.tsv

The inputs folder has example inputs example_InputGprime.tsv and example_InputGprime_1.tsv to execute and test the implementation.
