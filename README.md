# sc_longitudinal_infer

The Equation 2 is implemented in the file equation2Value.py and it accepts the following parameters as an input: the input data matrix, two output files from SCG (genotype_posterior.tsv and cluster_posterior.tsv), the lambda coefficient value and number of clusters. This file can be executed with the following command:

python equation2Value.py -input ../scg/examples/snv.tsv -gp ../no_doublet/genotype_posteriors.tsv -cp ../no_doublet/cluster_posteriors.tsv -lambda_coeff 0.2 -subclones 40

Please replace it with corresponding path, file names and parameters.
