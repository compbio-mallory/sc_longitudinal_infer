# sc_longitudinal_infer

For clustering we used SCG. Here are the steps to run SCG:

1. conda activate scg
2. scg run_singlet_model --config_file examples/config.yaml --out_dir ../no_doublet

The files generated are cluster_posteriors.tsv and genotype_posteriors.tsv.
The cluster numbers can be changed in the config.yaml file. 
---------------------------------------------------------------------------------------

The Equation 2 is implemented in the file equation2Value.py and it accepts the following parameters as an input: the input data matrix, two output files from SCG (genotype_posterior.tsv and cluster_posterior.tsv), the lambda coefficient value and number of clusters. This file can be executed with the following command:

python equation2Value.py -input ../scg/examples/snv.tsv -gp ../no_doublet/genotype_posteriors.tsv -cp ../no_doublet/cluster_posteriors.tsv -lambda_coeff 0.2 -subclones 40

Please replace it with corresponding path, file names and parameters.

---------------------------------------------------------------------------------

To run the longitudinal tree program use the following command: python longitudinalTree.py -input example_InputGprime.tsv

The inputs folder has example inputs example_InputGprime.tsv and example_InputGprime_1.tsv to execute and test the implementation.
