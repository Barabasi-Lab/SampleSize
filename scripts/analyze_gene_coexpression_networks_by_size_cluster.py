import sys, os
import configparser 
import hashlib
import optparse
import pwd
import subprocess


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_gene_coexpression_networks_by_size_cluster.py -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wto_N50 -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_sample_size -p 0.05
    """
    options = parse_options()
    analyze_gene_coexpression_networks_of_same_size(options)


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("analyze_gene_coexpression_networks_by_size_cluster.py  -s SIZE")

    # Directory arguments
    parser.add_option("-n", action="store", type="string", dest="networks_dir", help="Directory with the input networks", metavar="NETWORKS_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Directory for the output files and plots", metavar="OUTPUT_DIR")
    parser.add_option("-p", action="store", type="string", dest="pval_adj_cutoff", default=0.05, help="Cut-off of the adjusted p-value to consider an edge significant", metavar="PVAL_ADJ_CUTOFF")
    parser.add_option("-t", action="store", type="string", dest="top_percent", default=None, help="Top percentage of edges to select. If None, all of them are selected.", metavar="TOP_PERCENT")

    (options, args) = parser.parse_args()

    if options.networks_dir is None or options.output_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def analyze_gene_coexpression_networks_of_same_size(options):
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
    networks_dir = options.networks_dir
    output_dir = options.output_dir
    create_directory(output_dir)
    logs_dir = os.path.join(src_path, '../logs')
    create_directory(logs_dir)
    dummy_dir = os.path.join(src_path, '../cluster_scripts')
    create_directory(dummy_dir)
    script_file = os.path.join(src_path, 'analyze_networks_of_same_size.py')

    # Define queue parameters
    max_mem = config.get("Cluster", "max_mem")
    queue = config.get("Cluster", "cluster_queue") # debug, express, short, long, large
    modules = ['python/3.8.1', 'anaconda3']
    conda_environment = 'ScipherSampleSizeEnv'
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
    pval_adj_cutoff = float(options.pval_adj_cutoff)
    top_percent = options.top_percent
    repetitions = 10
    step = 10
    max_num_samples = 56
    #max_num_samples = 20
    size_list = list(range(10, max_num_samples+1, step))
    for size in size_list:

        if limit: # Break the loop if a limit of jobs is introduced
            if l > limit:
                print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                break

        bash_script_name = 'analysis_size_{}_pval_{}.sh'.format(size, pval_adj_cutoff)
        bash_script_file = os.path.join(dummy_dir, bash_script_name)

        output_file = os.path.join(output_dir, 'topology_size_{}.txt'.format(size))
        #if not fileExist(bash_script_file):
        if not fileExist(output_file):

            if top_percent:
                command = 'python {} -n {} -o {} -s {} -r {} -p {} -t {}'.format(script_file, networks_dir, output_dir, size, repetitions, pval_adj_cutoff, top_percent)
            else:
                command = 'python {} -n {} -o {} -s {} -r {} -p {}'.format(script_file, networks_dir, output_dir, size, repetitions, pval_adj_cutoff)
            print(command)
            submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name, constraint=True, conda_environment=conda_environment)

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

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, queue_parameters={'max_mem':5000, 'queue':'short', 'max_time':'24:00:00', 'logs_dir':'/tmp', 'modules':['Python/3.6.2']}, dummy_dir="/tmp", script_name=None, constraint=False, conda_environment=None):
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

