import sys, os
import configparser 
import hashlib
import optparse
import pwd
import subprocess


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_gene_coexpression_networks_cluster.py -i /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female -p /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female -m wgcna
    """
    options = parse_options()
    create_gene_coexpression_networks(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("analyze_gene_coexpression_networks_cluster.py")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Directory with the input datasets", metavar="NETWORKS_DIR")
    parser.add_option("-p", action="store", type="string", dest="ppi_file", help="File with the PPI network", metavar="PPI_FILE")
    parser.add_option("-d", action="store", type="string", dest="disease_genes_file", help="File with the disease-gene associations", metavar="DISEASE_GENES_FILE")
    parser.add_option("-e", action="store", type="string", dest="essential_genes_file", help="File with the essential genes", metavar="ESSENTIAL_GENES_FILE")
    parser.add_option("-g", action="store", type="string", dest="genes_dataset_file", help="File with the genes in the expression dataset", metavar="GENES_DATASET_FILE")
    parser.add_option("-o", action="store", type="string", dest="output_analysis_dir", help="Directory for the output analysis", metavar="OUTPUT_ANALYSIS_DIR")
    parser.add_option("-n", action="store", type="string", dest="output_networks_dir", help="Directory for the output networks", metavar="OUTPUT_NETWORKS_DIR")

    (options, args) = parser.parse_args()

    if options.input_dir is None or options.ppi_file is None or options.disease_genes_file is None or options.essential_genes_file is None or options.genes_dataset_file is None or options.output_analysis_dir is None or options.output_networks_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_gene_coexpression_networks(options):
    """
    Creates gene co-expression networks using the cluster
    """

    # Add "." to sys.path #
    src_path =  os.path.abspath(os.path.dirname(__file__))
    sys.path.append(src_path)

    # Read configuration file #     
    config = configparser.ConfigParser()
    config_file = os.path.join(src_path, "config.ini")
    config.read(config_file)

    # Define directories and files
    data_dir = os.path.join(src_path, '../data')
    input_dir = options.input_dir
    ppi_file = options.ppi_file
    disease_genes_file = options.disease_genes_file
    essential_genes_file = options.essential_genes_file
    genes_dataset_file = options.genes_dataset_file
    output_analysis_dir = options.output_analysis_dir
    create_directory(output_analysis_dir)
    output_networks_dir = options.output_networks_dir
    create_directory(output_networks_dir)
    logs_dir = os.path.join(src_path, '../logs')
    create_directory(logs_dir)
    dummy_dir = os.path.join(src_path, '../cluster_scripts')
    create_directory(dummy_dir)

    # Define queue parameters
    #max_mem = config.get("Cluster", "max_mem")
    max_mem = 60000
    queue = config.get("Cluster", "cluster_queue") # debug, express, short, long, large
    constraint = False
    exclude = []
    #exclude = ['d0012', 'd0123']
    #modules = ['R/4.2.0']
    modules = ['singularity/3.5.3']
    #conda_environment = 'SampleSizeR'
    conda_environment = None
    rstudio_environment = "/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"
    max_time_per_queue = {
        'debug'   : '0:20:00',
        'express' : '0:60:00',
        #'short'   : '24:00:00',
        'short'   : '5:00:00',
        'long'    : '5-0',
        'large'   : '6:00:00'
    }
    queue_parameters = {'max_mem':max_mem, 'queue':queue, 'max_time':max_time_per_queue[queue], 'logs_dir':logs_dir, 'modules':modules}

    # Limit of jobs
    limit = int(config.get("Cluster", "max_jobs_in_queue"))
    l = 1

    # Run co-expression for all files
    datasets = [f for f in os.listdir(input_dir) if fileExist(os.path.join(input_dir, f))]
    sizes = [str(size) for size in range(20, 12000, 20)]
    #sizes = [str(10)] + [str(size) for size in range(100, 12000, 100)]
    reps = [str(rep) for rep in range(1, 11, 1)]
    #reps = ['1', '2', '3']
    #reps = ['1']

    for dataset in sorted(datasets):

        dataset_short = dataset.replace('.net', '')
        dataset_info = dataset_short.split("_")
        method = dataset_info[0]
        size = dataset_info[-3]
        rep = dataset_info[-1]

        script_file = os.path.join(src_path, 'analyze_coexpression_network_by_significant_edges.R')
        coexpression_file = os.path.join(input_dir, dataset)
        
        dataset_names = []
        if size in sizes and rep in reps:
            if (method == "aracne"):
                threshold = 0
                dataset_name = '{}_threshold_{}'.format(dataset_short, str(threshold))
                dataset_names.append((dataset_name, threshold))
            else:
                #thresholds = [0.05, 0.001]
                thresholds = [0.05]
                for threshold in thresholds:
                    dataset_name = '{}_threshold_{}'.format(dataset_short, str(threshold))
                    dataset_names.append((dataset_name, threshold))
        
            for dataset_name, threshold in dataset_names:

                if limit: # Break the loop if a limit of jobs is introduced
                    if l > limit:
                        print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                        sys.exit(0)
                
                bash_script_name = 'analyze_coexpression_network_by_significant_edges_{}.sh'.format(dataset_name)
                bash_script_file = os.path.join(dummy_dir, bash_script_name)
                create_directory(os.path.join(output_networks_dir, 'subgraphs'))
                create_directory(os.path.join(output_networks_dir, 'main_core'))
                create_directory(os.path.join(output_networks_dir, 'ppi'))
                create_directory(os.path.join(output_networks_dir, 'ppi_main_core'))
                create_directory(os.path.join(output_networks_dir, 'disease_genes'))
                create_directory(os.path.join(output_networks_dir, 'essential_genes'))
                output_disease_genes_subgraph_dir = os.path.join(output_networks_dir, 'disease_genes/{}_{}'.format(dataset_name, 'disease_genes_subgraphs'))
                create_directory(output_disease_genes_subgraph_dir)
                output_topology_file = os.path.join(output_analysis_dir, '{}_{}'.format(dataset_name, 'analysis_topology.txt'))
                output_ppi_file = os.path.join(output_analysis_dir, '{}_{}'.format(dataset_name, 'analysis_ppi.txt'))
                output_disease_genes_file = os.path.join(output_analysis_dir, '{}_{}'.format(dataset_name, 'analysis_disease_genes.txt'))
                output_essential_genes_file = os.path.join(output_analysis_dir, '{}_{}'.format(dataset_name, 'analysis_essential_genes.txt'))

                if ((method in ["pearson", "spearman", "wto", "aracne"]) and (threshold in [0, 0.001, 0.05])):
                #if ((method in ["pearson", "spearman", "aracne", "wto"]) and (threshold in [0, 0.05])):
                #if ((method in ["pearson", "spearman"]) and (threshold in [0, 0.05])):

                    #if not fileExist(bash_script_file):
                    if ((not fileExist(output_topology_file)) or (not fileExist(output_ppi_file)) or (not fileExist(output_disease_genes_file)) or (not fileExist(output_essential_genes_file))):

                        command = 'Rscript {} -c {} -o {} -s {} -f {} -t {} -p {} -d {} -e {} -g {}'.format(script_file, coexpression_file, output_analysis_dir, output_networks_dir, dataset_name, threshold, ppi_file, disease_genes_file, essential_genes_file, genes_dataset_file)
                        print(dataset_name, threshold)
                        #print(output_topology_file)
                        print(command)
                        print(l)
                        submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=constraint, exclude=exclude, conda_environment=conda_environment, rstudio_environment=rstudio_environment)

                        l += 1

        # if limit: # Break the loop if a limit of jobs is introduced
        #     if l > limit:
        #         print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
        #         break

        # script_file = os.path.join(src_path, 'compare_coexpression_to_ppi.R')
        # coexpression_file = os.path.join(input_dir, dataset)
        # dataset_name = '{}'.format('.'.join(dataset.split('.')[:-1]))
        # output_file = os.path.join(output_analysis_dir, 'comparison_coexpression_ppi_{}.txt'.format(dataset_name))
        # bash_script_name = 'comparison_coexpression_ppi_{}.sh'.format(dataset_name)
        # bash_script_file = os.path.join(dummy_dir, bash_script_name)

        # #if not fileExist(bash_script_file):
        # if not fileExist(output_file):

        #     if 'aracne' not in dataset_name:
        #         # Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/compare_coexpression_to_ppi.R -c /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -p /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net/comparison_coexpression_ppi.txt
        #         command = 'Rscript {} -c {} -p {} -o {}'.format(script_file, coexpression_file, ppi_file, output_file)
        #         print(command)
        #         submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=constraint, exclude=exclude, conda_environment=conda_environment)

        #         l += 1

        # if limit: # Break the loop if a limit of jobs is introduced
        #     if l > limit:
        #         print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
        #         break
        
        # # Run the analysis of stability
        # script_file = os.path.join(src_path, 'analyze_stability_coexpression_networks.R')
        # coexpression_file = os.path.join(input_dir, dataset)
        # # Check if there is the word size in the name of the dataset (meaning it is not the network from all samples)
        # if 'size' in dataset:
        #     dataset_all_samples = dataset.split('size')[0] + 'all_samples.net'
        #     coexpression_file_all_samples = os.path.join(input_dir, dataset_all_samples)
        #     dataset_name = '{}'.format('.'.join(dataset.split('.')[:-1]))
        #     output_difference_file = os.path.join(output_analysis_dir, 'analysis_score_difference_{}.txt'.format(dataset_name))
        #     output_scores_file = os.path.join(output_analysis_dir, 'analysis_score_ranges_{}.txt'.format(dataset_name))
        #     output_threshold_file = os.path.join(output_analysis_dir, 'analysis_score_thresholds_{}.txt'.format(dataset_name))
        #     output_disease_file = os.path.join(output_analysis_dir, 'analysis_score_diseases_{}.txt'.format(dataset_name))
        #     bash_script_name = 'analysis_stability_coexpression_networks_{}.sh'.format(dataset_name)
        #     bash_script_file = os.path.join(dummy_dir, bash_script_name)

        #     # Check that the network of all samples exist!
        #     if fileExist(coexpression_file_all_samples):
        #         #if not fileExist(bash_script_file):
        #         if not fileExist(output_disease_file):

        #             if 'aracne' not in dataset_name:
        #                 # Rscript /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_stability_coexpression_networks.R -c /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.net -a /scratch/j.aguirreplans/Scipher/SampleSize/networks_GTEx/Spleen_female/wgcna_RNAseq_samples_Spleen_female_all_samples.net -p /home/j.aguirreplans/data/PPI/interactome_tissue_specific/interactome_2019_Spleen_female.csv -d /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/analysis_score_difference_wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.txt -s /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/analysis_score_ranges_wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.txt -t /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/networks_GTEx/Spleen_female/analysis_score_thresholds_wgcna_RNAseq_samples_Spleen_female_size_10_rep_2.txt
        #                 command = 'Rscript {} -c {} -a {} -p {} -d {} -s {} -t {} -r {}'.format(script_file, coexpression_file, coexpression_file_all_samples, ppi_file, output_difference_file, output_scores_file, output_threshold_file, output_disease_file)
        #                 print(command)
        #                 submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=constraint, exclude=exclude, conda_environment=conda_environment)

        #                 l += 1

    return




#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################


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

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, queue_parameters={'max_mem':5000, 'queue':'short', 'max_time':'24:00:00', 'logs_dir':'/tmp', 'modules':['Python/3.6.2']}, dummy_dir="/tmp", script_name=None, constraint=False, exclude=[], conda_environment=None, rstudio_environment=None):
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
                if rstudio_environment:
                    fd.write('RSTUDIO_IMAGE="{}"\n'.format(rstudio_environment))
                    fd.write('TMP=$HOME/tmp\n')
                    command = 'singularity run -B "/scratch:/scratch,/work:/work,${PWD}/tmp:/tmp" $RSTUDIO_IMAGE ' + command
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

