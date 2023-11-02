# scLongTree
Computational tool to infer longitudinal tree for scDNAseq data

## Table of Contents
- [Commands to run BnpC, scLongTree, getting SNV placements for all tools, evaluation and simulation.](#commands_3methods)
    * [BnpC](#bnpc)
    * [scLongTree](#scLongTree)
    * [SNV placements](#snvAccuracy)
    * [Evaluation](#evaluation)
    * [Simulation](#simulation)

# <a name="usage_of_scDNAseq_clustering"></a>Usage of scDNAseq Clustering.
## <a name="software_requirements"></a>Software Requirements ##

1. Python 2.7.15 or up.

# <a name="commands_3methods"></a>Commands to run BnpC, scLongTree, getting SNV placements for all tools, evaluation and simulation. #

## <a name="bnpc"></a>BnpC ##

Since we used BnpC as a clustering algorithm in our pipeline for the evaluation we will briefly list the steps we followed:

1. We process the noisy input matrix provided in ``` sample/input.D.csv ``` with the following scripts provided in ``` Clustering ```

	``` python processSimInput.py -input $j -output outputFile.csv ```

2. We run Bnpc on the processed input data and get the consensus genotype using the following script:

	``` python gen_bnpc_getGmatrix.py -cc assignment.txt -gp BnpCs_genotypes_posterior_mean.tsv -celltp sample/cellsInSubclones.csv -op consensus_genotype.tsv ```

3. You can also run your own clustering algorithm and get a cluster assignment file similar to assignment.txt and a consensus genotype similar to one provided in ``` Clustering```

## <a name="scLongTree"></a>scLongTree ## 

To run the algorithm we have all the necessary scripts in ``` algorithm ``` and sample inputs in ``` sample ```. The below steps will be explained with scripts provided in these folders.

1. Once we have the clustering results we will stratify it across timepoints with the following script:

	``` python algorithm/inputForTree.py -cc Clustering/assignment.txt -cg Clustering/consensus_genotype.tsv -sim false -celltp sample/cellsInSubclones.csv -op consensus_genotype_withTimepoints.csv -ccop cluster_cell.npy ```

2. Next we run our algorithm to infer the longitudinal tree:

	``` python algorithm/longitudinalTree.py -input consensus_genotype_withTimepoints.csv -cg Clustering/consensus_genotype.tsv -cc cluster_cell.npy -D sample/input.D.csv -celltp sample/cellsInSubclones.csv -sim false -op_tree inferredTree.csv -cloneNode cloneNode.npy -cloneCells cloneCells.npy ```

3. Now, given that we run the algorithm multiple times we select the best tree with following script. You can pass n no. of trees as an argument. I passed three trees below as arguments:

	``` python selectBestClusterAndTree.py -tree inferredTree1.csv inferredTree2.csv inferredTree3.csv ```

## <a name="snvAccuracy"></a>Calculating SNV accuracy ##

1. Given the results we evaluate based on placement of SNVs on the tree. Since the tree inferred is different for each tool we used the following scripts for each respective trees to get the pair of SNVS for ancestral, on same branch and incomparable/parallel:

scLongTree: ``` python algorithm/getSNVrelations.py -tree inferredTree.csv -output outputFileLoc ```	
LACE:	``` python LACE_exp/getSNVrelations.py -B lace_B.json -output outputFileLoc ``` 
SCITE: ``` python SCITE_exp/getSNVpairs.py -tree sciteTree.gv -output outputFileLoc ```
SiCloneFit: ``` python SiCloneFit_exp/getSNVpairs.py -tree best_MAP_tree.txt -genotype best_MAP_predicted_genotype.txt -output outputFileLoc ```

## <a name="evaluation"></a>Evaluation ##

1. Given the pair of SNVs we get the SNV pair accuracy by running the following script:

	``` python evaluation/evaluateSNVpairs.py -gt groundTruth_locationOf_SavedSNVpairs -scLongTree locationOf_SavedSNVpairs -siclonefit locationOf_SavedSNVpairs -scite locationOf_SavedSNVpairs -lace locationOf_SavedSNVpairs ```

2. We calculate the precision and recall for LACE and scLongTree:

	``` python calculate_PrecisionRecall.py -mut sample/clone_mut.csv -gtTree groundTruthTree -lgTree inferredTree -lace true/false ```
algorithm - has the recent implemented longitudinal algorithm

## <a name="simulation"></a>Simulation ##

### Steps to run the simulator using different variables. ###

1. To generate the tree with beta-splitting variable and desired number of leaves use the following script:
        * ``` python gen_tree.py -F $i -B 0.2 -o output_tree.csv ```
        * Parameters:
                * ```-F``` pass the number of leaves or desired number of clones.
                * ```-B``` pass the beta split variable value.
                * ```-o``` output tree name to save.

2. Once we have the tree we can assign cells, mutations and include false positives, false negatives, missing data, doublets using the following script:
        * ``` python sim_par.py -a 0.01 -b 0.2 -m 0.2 -c 500 -n 200 -e 0.1 -f output_tree.csv -P outputFilePrefix ```
        * Parameters:
                * ```-a``` false positive rate.
                * ```-b``` false negative rate.
                * ```-m``` missing rate.
                * ```-c``` number of cells.
                * ```-n``` number of mutations.
                * ```-e``` doublet rate.
                * ```-f``` the tree generated earlier using number of leaves and beta splitting variable.
                * ```-P``` the prefix and directory for output files.

Clustering - has the clustering pipeline with BnpC

evaluation - has the evaluation scripts to calculate precision, recall, V-measure, unobserved subclones.

inputs - has some old dummy data.

LACE_exp - has the pipeling to run LACE on our simulated data.

real_data_exp - has the pipeling to run scLongTree on real data SA501.

simulated_data_exp - has the pipeline to run scLongTree on our simulated data.

simulation - script to generate simulated data and all the data.

visualize - scripts to visualize the tree.

-------------------------------------------------------------------------------

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

UPDATE: The folders LACE_SIM_DATA* has the correct and updated equation2Value.py.

---------------------------------------------------------------------------------

Longitudinal Tree program accepts both matrix as tsv and a json as an input. To run the Longitudinal Tree program with a matrix as input use the following command:

python longitudinalTree.py -input ./examples/example_InputGprime.tsv -inputType matrix

To run with a dictionary(json) where subclones are already arranged in timepoints use the following command:

python longitudinalTree.py -input ./examples/simulated_input_tree.json -inputType dict

The inputs folder has example inputs example_InputGprime.tsv and example_InputGprime_1.tsv to execute and test the implementation.

-----------------------------------------------------------------------------------

To test with the simulated data or examples created by us run the scripts in this folder: our_simData_exp

------------------------------------------------------------------------------------

<u> Steps to run experiments on LACE Breast cancer real data: <u> 

1. Scripts to prepare the data for testing from Rdata to json or CSV can be found in /gpfs/research/fangroup/rk18g/longitudinal/ and steps are explained in the processData.txt (two files to run: convertRToCSV.R and convertRToJson.R)
2. Check any of the folders starting with lace_data_exp* inside /gpfs/research/fangroup/rk18g/longitudinal/SCG_Roth/Gyanendra_code/
3. We have the input files as merged.tsv and should run the following set of commands in sequence from inside the folder lace_data_exp*:
  
   python save_multipleSCGresults.py -input merged.tsv.gz -cluster_list 5 6 7 8 9 10 -scg_config ../../../SCG_Roth/scg/examples/config.yaml -niters 10 (Run this file from the SCG conda environment)
  
  python chooseMinVal_Eq2.py -input merged.tsv.gz -lambda_coeff 0.2

  Next get the timepoints for each cell using the following script:
/gpfs/research/fangroup/rk18g/longitudinal/prepare_dataForClustering.py

  python input_for_tree.py -input merged.tsv.gz -gp ./scg_clusterNo_5_iter_6/genotype_posteriors.tsv.gz -cellCluster ./scg_clusterNo_5_iter_6/cellToCluster.json -timepoints lace_real_data_tp_1.json -output inputForLongitudinalTree.json

  python longitudinalTree.py -input inputForLongitudinalTree.json -inputType dict
  
  -------------------------------------------------------------------------------
  
 Steps to get LACE’s simulated data ready for further analysis:
  
1. For each topologies we have a file extractData.py which helps to extract each datasets from the topologies. Need to manually mention inside code which dataset we want to extract from the input_data.json.

2. For each extracted dataset in json formats we execute a file processDatasets.py which helps to get the raw D matrix, ground truth matrix and their provided corrected genotype matrix. Modify the # of cells in each time point according to the dataset.

3. After this we use the input_data.rdata and execute the extractData_usingR.R in RStudio to get lace_inference.rdata. Once we have that we use convertRToJson.R to convert the inference results to a json file.

4. Then we use analyze_inference.py script to get the LACE_inference results as matrix. Here also we have to change # of cells in each timepoint. 

5. Copy the D, ground truth and lace inference matrix to respective folders for downstream analysis.
  
----------------------------------------------------------------------------------
  
  Steps to run experiments with LACE Simulated data:

1. Inside /gpfs/research/fangroup/rk18g/longitudinal/LACE-UTILITIES/simulations/exp1_smartseq_10x/topology_1 there are two scripts: extractData.py and processDatasets.py. The first one gets sample dataset from each topology and second script helps to get the input data D (DMatrix_1.tsv), ground truth matrix (GT_matrix_1.tsv) and corrected genotype matrix.
2. We can copy the required matrices in the following path: /gpfs/research/fangroup/rk18g/longitudinal/SCG_Roth/Gyanendra_code/lace_sim_data.
  Then follow the following steps:
  
  1. python save_multipleSCGresults.py -input DMatrix_1.tsv.gz -cluster_list 20 40 60 80 100 -scg_config ../../../SCG_Roth/scg/examples/config.yaml -niters 10 (to run SCG with many iterations).
  2. python chooseMinVal_Eq2.py -input DMatrix_1.tsv.gz -lambda_coeff 0.2 (choose the SCG result with min Eq 2 value). (This also returns the cell timepoints)
  3. python input_for_tree.py -input DMatrix_1.tsv.gz -gp ./scg_clusterNo_20_iter_8/genotype_posteriors.tsv.gz -cellCluster ./scg_clusterNo_20_iter_8/cellToCluster.json -timepoints cell_timepoints.json -output inputForLongitudinalTree.json
  4. python longitudinalTree.py -input inputForLongitudinalTree.json -inputType dict
  
  This will give you the tree using our algorithm. To compare with LACE's tree that is saved in MacBook's Documents: extractData_usingR.R and can be executed in RStudio.
  
3. The evaluation metrics which is to get the FP, TP, FN, TN values comapring LACE's data, inferred data and ground truth data can be done in the following way:
  1. python equation2Value.py -input DMatrix_1.tsv.gz -gp ./scg_clusterNo_20_iter_8/genotype_posteriors.tsv.gz -cp ./scg_clusterNo_20_iter_8/cluster_posteriors.tsv.gz -lambda_coeff 0.2 -subclones 20 -cellCluster ./scg_clusterNo_20_iter_8/cellToCluster.json (To get difference between D and our inG'', our inG'' and gtG'').
  2. python lace_equation2Value.py -inputD DMatrix_1.tsv.gz -inputG Lace_inferred_G.tsv  (Get LACE's results between inG'' and gtG'', D and inG'', D and gtG'')
  3. python evaluateMetrics.py -inputD DMatrix_1.tsv.gz -inG GDoublePrimeDf.tsv -gtG GT_matrix_1.tsv (Get corrected FP,FN, lost TP,TN, inferred FP,FN for both LACE and our data).

