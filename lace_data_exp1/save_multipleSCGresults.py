import yaml
import subprocess
import argparse

def read_yaml(yaml_file):
    with open(yaml_file,'r') as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        config_params = yaml.load(file, Loader=yaml.FullLoader)
        return config_params
        #print(config_params['num_clusters'])
        
def save_multipleSCGresults(config_params, fileName, clusterNo_list, nIters):
    for i in clusterNo_list:
        cluster_update = {'num_clusters': i}
        ufileName = {'data': {'snv': {'file': fileName, 'gamma_prior':[[98, 1, 1], [25, 50, 25], [1, 1, 98]], 'state_prior': [1, 1, 1]}}}
        #gamma_prior = {'data': {'snv': {'gamma_prior':[[98, 1, 1], [25, 50, 25], [1, 1, 98]]}}}
        #state_prior = {'data': {'snv': {'state_prior': [1, 1, 1]}}}
        config_params.update(cluster_update)
        config_params.update(ufileName)
        print(config_params)
        # Update the config file to have the values
        with open('../../SCG_Roth/scg/examples/new_config.yaml', 'w') as file:
            documents = yaml.dump(config_params, file)

        for j in range(1, nIters):
            print("Iteration ",j)
            dir_name = 'scg_clusterNo_'+str(i)+'_iter_'+str(j)
            subprocess.call('mkdir '+dir_name, shell=True)
            scg_cmd = 'scg run_singlet_model --config_file ../../SCG_Roth/scg/examples/new_config.yaml --out_dir '+dir_name
            subprocess.call(scg_cmd, shell=True)

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (similar to input to SCG)")
parser.add_argument("-cluster_list", "--cluster_list",dest = "cluster_list", help="List of cluster numbers")
parser.add_argument("-scg_config","--scg_config",dest="scg_config", help="SCG Config yaml file")
parser.add_argument("-niters", "--niters",dest ="niters", help="No of iterations")
args = parser.parse_args()

config_params = read_yaml(args.scg_config)
clusterNo_list = args.cluster_list
save_multipleSCGresults(config_params, args.input, clusterNo_list, args.niters)

#config_params = read_yaml('../../SCG_Roth/scg/examples/config.yaml')
#clusterNo_list = [5,6,7,8,9,10]
#save_multipleSCGresults(config_params, '/gpfs/research/fangroup/rk18g/longitudinal/LACE-UTILITIES/real_data/Breast_Cancer/p494/final_LACE_input_data/merged.tsv.gz',clusterNo_list, 10)

# update dictionary config_params and dump to the same yaml file
