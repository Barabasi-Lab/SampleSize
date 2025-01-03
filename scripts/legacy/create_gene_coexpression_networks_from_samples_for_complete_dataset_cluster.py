import sys, os
import configparser 
import hashlib
import optparse
import pwd
import subprocess


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_networks_from_samples_for_complete_dataset_cluster.py -i /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/GTEx/sampling_with_repetition -r /work/ccnr/j.aguirreplans/Databases/GTEx/v8/tpm_filtered_files_by_tissue -o /work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_gtex -m pearson
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_networks_from_samples_for_complete_dataset_cluster.py -i /home/j.aguirreplans/Projects/Scipher/SampleSize/data/sampling/TCGA/2022-07-27-Dataset/sampling_with_repetition/tumor -r /work/ccnr/j.aguirreplans/Databases/TCGA/2022-07-27-Dataset/TCGA/out/tumor/filter_genes_low_counts/rnaseq_filtered_files_by_project -o /work/ccnr/j.aguirreplans/Scipher/SampleSize/networks_tcga -m pearson
    """
    options = parse_options()
    create_gene_coexpression_networks(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("create_gene_coexpression_networks_cluster.py  -m METRIC")

    # Directory arguments
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Dataset directory with folders containing input sample files", metavar="NETWORKS_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Dataset directory to store the output networks", metavar="OUTPUT_DIR")
    parser.add_option("-r", action="store", type="string", dest="rnaseq_data", help="Dataset directory containing RNA-seq expression files. If -u included, then it's a RNA-seq file instead of a directory.", metavar="RNASEQ_DATA")
    parser.add_option("-m", action="store", type="string", dest="metric", default="pearson", help="Gene co-expression method", metavar="METRIC")
    parser.add_option("-u", action="store_true", dest="unique_rnaseq_file", help="If included, the option -r is a file that will be used for all the datasets")

    (options, args) = parser.parse_args()

    if options.input_dir is None or options.output_dir is None or options.rnaseq_data is None or options.metric is None:
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

    metric = options.metric
    wto_n = 100
    wto_delta = 0.2
    wgcna_power = 6
    wgcna_type = "signed"
    mi_estimator = "pearson"
    aracne_eps = 0
    correction_method = "bonferroni"

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
    rnaseq_data = options.rnaseq_data
    output_dir = options.output_dir
    create_directory(output_dir)
    logs_dir = os.path.join(src_path, '../logs')
    create_directory(logs_dir)
    dummy_dir = os.path.join(src_path, '../cluster_scripts')
    create_directory(dummy_dir)
    #script_file = os.path.join(src_path, 'test.R')
    script_file = os.path.join(src_path, 'create_gene_coexpression_network_from_samples_list.R')

    # Define queue parameters
    max_mem = config.get("Cluster", "max_mem")
    queue = config.get("Cluster", "cluster_queue") # debug, express, short, long, large
    constraint = True
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

    # Run co-expression for all types of datasets
    #dataset_types_selected = ["TCGA-BRCA", "TCGA-LGG", "TCGA-LUSC", "TCGA-LUAD", "TCGA-THCA", "TCGA-COAD", "TCGA-PRAD", "Breast.Mammary.Tissue", "Brain.Cortex", "Lung", "Thyroid", "Colon.Transverse", "Colon.Sigmoid", "Prostate"]
    #dataset_types_selected = ["TCGA-BRCA", "TCGA-THCA", "TCGA-LGG", "TCGA-LUSC", "TCGA-COAD"]
    dataset_types_selected = ["Breast.Mammary.Tissue", "Thyroid", "Lung", "Colon.Transverse", "Colon.Sigmoid"]
    #dataset_types = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    dataset_types = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d)) and d in dataset_types_selected]
    #sizes = [str(size) for size in range(10, 12000, 20)]
    #sizes = [str(size) for size in range(20, 12000, 20)]
    sizes = [10] + [str(size) for size in range(20, 12000, 20)]
    #sizes = [10] + [str(size) for size in range(100, 12000, 100)]
    #reps = [str(rep) for rep in range(5, 11, 1)]
    #reps = [str(rep) for rep in range(1, 11, 1)]
    reps = ['1', '2', '3', '4', '5']
    #reps = ['1', '2']

    for dataset_type in sorted(dataset_types):

        datasets = sorted([f for f in os.listdir(os.path.join(input_dir, dataset_type)) if fileExist(os.path.join(input_dir, dataset_type, f))])

        # Get rnaseq file
        if options.unique_rnaseq_file:
            rnaseq_file = rnaseq_data
        else:
            test_file = [f for f in os.listdir(rnaseq_data) if fileExist(os.path.join(rnaseq_data, f))][0]
            prefix_file = '_'.join(test_file.split('_')[:-1])
            sufix_file = test_file.split('.')[-1]
            rnaseq_file = os.path.join(rnaseq_data, '{}_{}.{}'.format(prefix_file, dataset_type, sufix_file))
        
        # Get output dir
        dataset_output_dir = os.path.join(output_dir, dataset_type)
        create_directory(dataset_output_dir)

        for x in range(len(datasets)):

            dataset = datasets[x]
            size = dataset.replace('.txt', '').split("_")[-3]
            rep = dataset.replace('.txt', '').split("_")[-1]

            # Check if RNAseq exists
            if not fileExist(rnaseq_file):
                print('RNAseq file "{}" does not exist for the corresponding dataset "{}"'.format(rnaseq_file, dataset))
                sys.exit(10)

            if ((size in sizes) and (rep in reps)):
            #if ((size in sizes) and (rep in reps) and (size not in ['60', '80'])):

                if limit: # Break the loop if a limit of jobs is introduced
                    if l > limit:
                        print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                        break

                samples_file = os.path.join(input_dir, dataset_type, dataset)
                dataset_name = '{}_{}'.format(metric, '.'.join(dataset.split('.')[:-1]))
                output_file = os.path.join(dataset_output_dir, '{}.net'.format(dataset_name))
                bash_script_name = '{}.sh'.format(dataset_name)
                bash_script_file = os.path.join(dummy_dir, bash_script_name)

                #if not fileExist(bash_script_file):
                if not fileExist(output_file):

                    #command = 'Rscript {}'.format(script_file)
                    command = 'Rscript {} -s {} -f {} -o {} -m {} -n {} -d {} -p {} -t {} -e {} -a {} -c {}'.format(script_file, samples_file, rnaseq_file, output_file, metric, wto_n, wto_delta, wgcna_power, wgcna_type, mi_estimator, aracne_eps, correction_method)
                    print(command)
                    print(size, rep)
                    print(l)
                    submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=constraint, exclude=exclude, conda_environment=conda_environment)

                    l += 1

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

