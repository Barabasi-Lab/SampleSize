import sys, os
import configparser 
import hashlib
import itertools
import optparse
import pandas as pd
import pwd
import subprocess


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_differentially_coexpressed_genes_cluster.py -d /work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_tcga/tumor/TCGA-BRCA -n /work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_gtex/Breast.Mammary.Tissue -o /work/ccnr/j.aguirreplans/Scipher/SampleSize/differential_coexpression_analysis/TCGA-BRCA_Breast.Mammary.Tissue -p 0.05
    """
    options = parse_options()
    calculate_differentially_coexpressed_genes(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("calculate_differentially_coexpressed_genes_cluster.py  -s SIZE")

    # Directory arguments
    parser.add_option("-d", action="store", type="string", dest="networks_dir_D", help="Directory with the input networks of condition 1 (e.g., disease)", metavar="NETWORKS_DIR_D")
    parser.add_option("-n", action="store", type="string", dest="networks_dir_N", help="Directory with the input networks of condition 2 (e.g., normal)", metavar="NETWORKS_DIR_N")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Directory for the output files", metavar="OUTPUT_DIR")
    parser.add_option("-p", action="store", type="float", dest="pval_adj_cutoff", default=0.05, help="Cut-off of the adjusted p-value to consider an edge significant", metavar="PVAL_ADJ_CUTOFF")
    parser.add_option("-s", action="store_true", dest="stretch_normalization", help="If included, uses stretch parameter to normalize networks before running CODINA")
    parser.add_option("-f", action="store_true", dest="filter_by_common_nodes", help="If included, filters networks by keeping nodes that are common in both networks")

    (options, args) = parser.parse_args()

    if options.networks_dir_D is None or options.networks_dir_N is None or options.output_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def calculate_differentially_coexpressed_genes(options):
    """
    Sends jobs that calculate differentially co-expressed genes between pairs of networks
    """

    # Add "." to sys.path #
    src_path =  os.path.abspath(os.path.dirname(__file__))
    sys.path.append(src_path)

    # Read configuration file #     
    config = configparser.ConfigParser()
    config_file = os.path.join(src_path, "config.ini")
    config.read(config_file)

    # Define directories and files
    networks_dir_D = options.networks_dir_D
    networks_dir_N = options.networks_dir_N
    output_dir = options.output_dir
    create_directory(output_dir)
    logs_dir = os.path.join(src_path, '../logs')
    create_directory(logs_dir)
    dummy_dir = os.path.join(src_path, '../cluster_scripts')
    create_directory(dummy_dir)
    script_file = os.path.join(src_path, 'calculate_differentially_coexpressed_genes.R')

    # Define queue parameters
    max_mem = config.get("Cluster", "max_mem")
    queue = config.get("Cluster", "cluster_queue") # debug, express, short, long, large
    constraint = False
    exclude = []
    #exclude = ['d0012', 'd0123']
    modules = ['python/3.8.1', 'anaconda3', 'R/4.0.3']
    #conda_environment = 'SampleSizeR'
    conda_environment = None
    max_time_per_queue = {
        'debug'   : '0:20:00',
        'express' : '0:60:00',
        'short'   : '24:00:00',
        'long'    : '5-0',
        'large'   : '6:00:00'
    }
    queue_parameters = {'max_mem':max_mem, 'queue':queue, 'max_time':max_time_per_queue[queue], 'logs_dir':logs_dir, 'modules':modules}

    # Limit of jobs
    limit = int(config.get("Cluster", "max_jobs_in_queue"))
    l = 1

    # Define additional parameters
    if options.stretch_normalization:
        stretch_normalization = '-s'
    else:
        stretch_normalization = ''
    
    if options.filter_by_common_nodes:
        filter_by_common_nodes = '-f'
    else:
        filter_by_common_nodes = ''

    # Get all possible files in both directories
    pval_adj_cutoff = float(options.pval_adj_cutoff)
    networks_D = [f for f in os.listdir(networks_dir_D) if fileExist(os.path.join(networks_dir_D, f))]
    networks_N = [f for f in os.listdir(networks_dir_N) if fileExist(os.path.join(networks_dir_N, f))]
    networks_df = pd.concat([get_info_from_network_files(networks_D, dataset_name="D"), get_info_from_network_files(networks_N, dataset_name="N")])
    networks_df['size'] = networks_df['size'].astype(int)
    sizes_d = set(networks_df[networks_df['dataset_name']=='D']['size'].tolist())
    sizes_n = set(networks_df[networks_df['dataset_name']=='N']['size'].tolist())
    shared_sizes = sorted(sizes_d & sizes_n)
    reps_d = set(networks_df[networks_df['dataset_name']=='D']['rep'].tolist())
    reps_n = set(networks_df[networks_df['dataset_name']=='N']['rep'].tolist())
    shared_reps = sorted(reps_d & reps_n)
    all_combinations = []
    for size in shared_sizes:
        # Select any combinations between networks of same size
        networks_d_selected = networks_df[(networks_df['size']==size) & (networks_df['dataset_name']=='D')]['network'].tolist()
        networks_n_selected = networks_df[(networks_df['size']==size) & (networks_df['dataset_name']=='N')]['network'].tolist()
        unique_combinations = list(itertools.product(networks_d_selected, networks_n_selected))
        all_combinations = all_combinations + unique_combinations

        # Select networks that have both same size and repetition number
        #for rep in shared_reps:
        #    network_d_selected = networks_df[(networks_df['size']==size) & (networks_df['rep']==rep) & (networks_df['dataset_name']=='D')]['network'].tolist()
        #    network_n_selected = networks_df[(networks_df['size']==size) & (networks_df['rep']==rep) & (networks_df['dataset_name']=='N')]['network'].tolist()
        #    if (len(network_d_selected) == 1) & (len(network_n_selected) == 1):
        #        #print([network_d_selected[0], network_n_selected[0]])
        #        all_combinations = all_combinations + [[network_d_selected[0], network_n_selected[0]]]

    for combination in all_combinations:

        if limit: # Break the loop if a limit of jobs is introduced
            if l > limit:
                print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                break

        network_name_d = combination[0]
        network_name_n = combination[1]
        bash_script_name = 'diffanalysis_{}___{}_pval_{}.sh'.format(network_name_d, network_name_n, pval_adj_cutoff)
        bash_script_file = os.path.join(dummy_dir, bash_script_name)

        output_edges_file = os.path.join(output_dir, 'diffanalysis_edges_{}___{}_pval_{}.txt'.format(network_name_d, network_name_n, pval_adj_cutoff))
        output_nodes_file = os.path.join(output_dir, 'diffanalysis_nodes_{}___{}_pval_{}.txt'.format(network_name_d, network_name_n, pval_adj_cutoff))
        #if not fileExist(bash_script_file):
        if not fileExist(output_edges_file) or not fileExist(output_nodes_file):

            command = 'Rscript {} -d {} -n {} -l {} -v {} -t {} {} {}'.format(script_file, os.path.join(networks_dir_D, network_name_d), os.path.join(networks_dir_N, network_name_n), output_edges_file, output_nodes_file, pval_adj_cutoff, stretch_normalization, filter_by_common_nodes)
            print(command)
            print(l)
            submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=constraint, exclude=exclude, conda_environment=conda_environment)

            l += 1

    return




#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


def get_info_from_network_files(network_files_list, dataset_name="D"):
    """
    Get the information about the size and repetitions from a list of networks
    """
    # Create a data frame to store the results
    networks_df = pd.DataFrame(columns=["network", "dataset_name", "size", "rep"])
    for network in network_files_list:
        network_name = '.'.join(network.split('.')[:-1])
        elem = network_name.split('_')
        if (elem[-1] == "consensus"):
            size = elem[-2]
            rep = "consensus"
        else:
            size = elem[-3]
            rep = elem[-1]
        # Add the information to the main data frame
        networks_df = networks_df.append(pd.DataFrame([[network, dataset_name, size, rep]], columns=["network", "dataset_name", "size", "rep"], index=[network_name]))
    return networks_df

def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return

#-------------#
# Cluster     #
#-------------#

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, queue_parameters={'max_mem':5000, 'queue':'short', 'max_time':'24:00:00', 'logs_dir':'/tmp', 'modules':['Python/3.6.2']}, dummy_dir="/tmp", script_name=None, constraint=False, exclude=[], conda_environment=None):
    """
    This function submits any {command} to a cluster {queue}.

    @input:
    command {string}
    queue {string} by default it submits to any queue
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue
    constraint is a boolean to say if you want to calculate the computations on specific nodes

    """
    #if max_jobs_in_queue is not None:
    #    while number_of_jobs_in_queue() >= max_jobs_in_queue: time.sleep(5)

    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    if not script_name:
        script_name = "submit_"+hashlib.sha224(command).hexdigest()+".sh"
    script= os.path.join(dummy_dir, script_name)
    if queue_file is not None:
        # Use a queue file as header of the bash file
        fd=open(script,"w")
        with open(queue_file,"r") as queue_standard:
            data=queue_standard.read()
            fd.write(data)
            fd.write("%s\n\n"%(command))
        fd.close()
        queue_standard.close()
        if queue is not None:
            os.system("sbatch -p %s %s" % (queue,script))
        else:
            os.system("sbatch %s"% (script))
    else:
        if queue_parameters is not None:
            # Write manually the bash file
            with open(script, "w") as fd:
                fd.write('#!/bin/bash\n') # bash file header
                fd.write('#SBATCH -p {}\n'.format(queue_parameters['queue'])) # queue name
                fd.write('#SBATCH --time={}\n'.format(queue_parameters['max_time'])) # max time in queue
                fd.write('#SBATCH --mem={}\n'.format(queue_parameters['max_mem'])) # max memory
                if constraint:
                    fd.write('#SBATCH --constraint="[cascadelake|zen2]"\n') # constraint to calculate on specific nodes
                if len(exclude) > 0:
                    fd.write('#SBATCH --exclude={}\n'.format(','.join(exclude))) # exclude specific nodes
                fd.write('#SBATCH -o {}.out\n'.format(os.path.join(queue_parameters['logs_dir'], '%j_{}'.format(script_name)))) # standard output
                fd.write('#SBATCH -e {}.err\n'.format(os.path.join(queue_parameters['logs_dir'], '%j_{}'.format(script_name)))) # standard error
                for module in queue_parameters['modules']:
                    fd.write('module load {}\n'.format(module)) # modules to load
                if conda_environment:
                    fd.write('source activate {}\n'.format(conda_environment))
                fd.write('{}\n'.format(command)) # command
            os.system("sbatch {}".format(script)) # execute bash file
        else:
            # Execute directly the command without bash file
            if queue is not None:
                os.system("echo \"%s\" | sbatch -p %s" % (command, queue))
            else:
                os.system("echo \"%s\" | sbatch" % command)


def number_of_jobs_in_queue():
    """
    This functions returns the number of jobs in queue for a given
    user.

    """

    # Initialize #
    user_name = get_username()

    process = subprocess.check_output(["squeue", "-u", user_name])

    return len([line for line in process.split("\n") if user_name in line])


def get_username():
    """
    This functions returns the user name.

    """

    return pwd.getpwuid(os.getuid())[0]


if  __name__ == "__main__":
    main()

