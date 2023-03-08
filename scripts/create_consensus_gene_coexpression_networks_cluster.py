import sys, os
import configparser 
import hashlib
import optparse
import pwd
import subprocess


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_consensus_gene_coexpression_networks_cluster.py
    """
    options = parse_options()
    create_consensus_gene_coexpression_networks_cluster(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("create_consensus_gene_coexpression_networks_cluster.py -i NETWORKS_DIR -o OUTPUT_DIR")

    # Directory arguments
    parser.add_option("-n", action="store", type="string", dest="networks_dir", help="Directory with the input networks", metavar="NETWORKS_DIR")
    #parser.add_option("-o", action="store", type="string", dest="output_dir", help="Directory for the output networks", metavar="OUTPUT_DIR")
    parser.add_option("-m", action="store", type="string", dest="method", default="pearson", help="File with the names of the drugs to analyze", metavar="METHOD")
    parser.add_option("-t", action="store", type="float", dest="threshold", default=0.05, help="File with the names of the drugs to analyze", metavar="THRESHOLD")

    (options, args) = parser.parse_args()

    if options.networks_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def create_consensus_gene_coexpression_networks_cluster(options):
    """
    Creates consensus gene co-expression networks using the cluster
    """

    # Add "." to sys.path #
    src_path =  os.path.abspath(os.path.dirname(__file__))
    sys.path.append(src_path)

    # Read configuration file #     
    config = configparser.ConfigParser()
    config_file = os.path.join(src_path, "config.ini")
    config.read(config_file)

    # Define directories and files
    networks_dir = options.networks_dir
    output_dir = os.path.join(networks_dir, 'consensus')
    create_directory(output_dir)
    inputs_dir = os.path.join(src_path, '../inputs')
    create_directory(inputs_dir)
    logs_dir = os.path.join(src_path, '../logs')
    create_directory(logs_dir)
    dummy_dir = os.path.join(src_path, '../cluster_scripts')
    create_directory(dummy_dir)
    script_file = os.path.join(src_path, 'create_consensus_gene_coexpression_network.R')

    # Define other inputs
    threshold = options.threshold
    method = options.method

    # Define queue parameters
    #max_mem = config.get("Cluster", "max_mem")
    max_mem = 120000
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
        'short'   : '24:00:00',
        'long'    : '5-0',
        'large'   : '6:00:00'
    }
    queue_parameters = {'max_mem':max_mem, 'queue':queue, 'max_time':max_time_per_queue[queue], 'logs_dir':logs_dir, 'modules':modules}

    # Limit of jobs
    limit = int(config.get("Cluster", "max_jobs_in_queue"))
    l = 1

    # Run co-expression for all files
    networks = [f for f in os.listdir(networks_dir) if fileExist(os.path.join(networks_dir, f))]
    #sizes = [str(size) for size in range(20, 12000, 20)]
    #sizes = [10] + [size for size in range(20, 12000, 20)]
    #sizes = [10] + [str(size) for size in range(100, 12000, 100)]
    #reps = [str(rep) for rep in range(5, 11, 1)]
    #reps = [str(rep) for rep in range(1, 11, 1)]
    #reps = ['1', '2']

    # Get network attributes
    size_to_networks = {}
    size_to_network_name = {}
    for network in sorted(networks):
        size = int(network.replace('.txt', '').split("_")[-3])
        network_name_without_rep = '_'.join(network.replace('.txt', '').split("_")[:-2])
        size_to_networks.setdefault(size, []).append(network)
        size_to_network_name[size] = network_name_without_rep
    
    for size in sorted(size_to_network_name.keys()):

        #if size in sizes:

            if limit: # Break the loop if a limit of jobs is introduced
                if l > limit:
                    print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                    break

            # Write list file
            network_name_without_rep = size_to_network_name[size]
            list_file = os.path.join(inputs_dir, 'network_replicates_{}.txt'.format(network_name_without_rep))
            with open(list_file, 'w') as list_fd:
                for network in size_to_networks[size]:
                    list_fd.write('{}\n'.format(os.path.join(networks_dir, network)))
            
            # Define other parameters
            output_file = os.path.join(output_dir, '{}_consensus.net'.format(network_name_without_rep))
            bash_script_name = '{}.sh'.format(network_name_without_rep)
            bash_script_file = os.path.join(dummy_dir, bash_script_name)

            #if not fileExist(bash_script_file):
            if not fileExist(output_file):

                #command = 'Rscript {}'.format(script_file)
                command = 'Rscript {} -l {} -n {} -t {} -m {}'.format(script_file, list_file, output_file, threshold, method)
                print(command)
                print(size)
                print(l)
                submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=constraint, exclude=exclude, conda_environment=conda_environment, rstudio_environment=rstudio_environment)

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

