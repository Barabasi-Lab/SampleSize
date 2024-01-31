import sys, os
import configparser 
import hashlib
import numpy as np
import pandas as pd
import pwd
import subprocess


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/create_gene_coexpression_networks_wto_cluster.py
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
    input_dir = os.path.join(data_dir, 'bootstrap/nonresponders')
    output_dir = os.path.join(data_dir, 'bootstrap_wTO/nonresponders')
    create_directory(output_dir)
    logs_dir = os.path.join(src_path, '../logs')
    create_directory(logs_dir)
    dummy_dir = os.path.join(src_path, '../cluster_scripts')
    create_directory(dummy_dir)
    script_file = os.path.join(src_path, '01_Generate_CN_with_wTO.R')

    # Define queue parameters
    max_mem = config.get("Cluster", "max_mem")
    queue = config.get("Cluster", "cluster_queue") # debug, express, short, long, large
    modules = ['R/4.0.3']
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
    datasets = [f for f in os.listdir(input_dir) if fileExist(os.path.join(input_dir, f))]

    for dataset in sorted(datasets):

        if limit: # Break the loop if a limit of jobs is introduced
            if l > limit:
                print('The number of submitted jobs arrived to the limit of {}. The script will stop sending submissions!'.format(limit))
                break

        input_file = os.path.join(input_dir, dataset)
        dataset_name = '.'.join(dataset.split('.')[:-1])
        output_file = os.path.join(output_dir, '{}.net'.format(dataset_name))
        bash_script_name = '{}.sh'.format(dataset_name)
        bash_script_file = os.path.join(dummy_dir, bash_script_name)

        #if not fileExist(bash_script_file):
        if not fileExist(output_file):

            command = 'Rscript {} -f {} -o {}'.format(script_file, input_file, output_file)
            print(command)
            submit_command_to_queue(command, max_jobs_in_queue=int(config.get("Cluster", "max_jobs_in_queue")), queue_file=None, queue_parameters=queue_parameters, dummy_dir=dummy_dir, script_name=bash_script_name)

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

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, queue_parameters={'max_mem':5000, 'queue':'short', 'max_time':'24:00:00', 'logs_dir':'/tmp', 'modules':['Python/3.6.2']}, dummy_dir="/tmp", script_name=None):
    """
    This function submits any {command} to a cluster {queue}.

    @input:
    command {string}
    queue {string} by default it submits to any queue
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue

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
                fd.write('#SBATCH -o {}.out\n'.format(os.path.join(queue_parameters['logs_dir'], script_name))) # standard output
                fd.write('#SBATCH -e {}.err\n'.format(os.path.join(queue_parameters['logs_dir'], script_name))) # standard error
                for module in queue_parameters['modules']:
                    fd.write('module load {}\n'.format(module)) # modules to load
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

