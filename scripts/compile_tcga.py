import sys, os
import json
import gzip
import numpy as np
import optparse
import pandas as pd


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/compile_tcga.py 
    """
    options = parse_options()
    compile_tcga(options)
    return


def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("analyze_gene_coexpression_networks_cluster.py")

    # Directory arguments
    parser.add_option("-d", action="store", type="string", dest="data_dir", default="/home/j.aguirreplans/Databases/TCGA/2022-03-28-Dataset/TCGA/raw/data", help="Directory with the input datasets", metavar="DATA_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_file", default="/home/j.aguirreplans/Databases/TCGA/2022-03-28-Dataset/TCGA/out/TCGA.csv", help="Path to the output file", metavar="OUTPUT_FILE")

    (options, args) = parser.parse_args()

    if options.data_dir is None or options.output_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#################
# MAIN FUNCTION #
#################

def compile_tcga(options):
    """
    Compile TCGA scattered datasets into a unique file
    """

    # Define main directories
    data_dir = options.data_dir
    output_file = options.output_file
    folder_list = os.listdir(data_dir)

    i=1
    if not fileExist(output_file):
        # Parse patients data
        main_data = None
        for patient_folder in folder_list:
            patient_folder_path = os.path.join(data_dir, os.fsdecode(patient_folder))
            if(os.path.isdir(patient_folder_path)):
                for file_id in os.listdir(patient_folder_path):
                    file_path = os.path.join(patient_folder_path, os.fsdecode(file_id))
                    if(file_id.split(".")[-1] == "gz" or file_id.split(".")[-1] == "counts"):
                        # Uncompress (if necessary) and read file
                        if(file_id.split(".")[-1] == "gz"):
                            f = gzip.open(file_path)
                            aux_data = pd.read_csv(f,sep='\t',index_col=0,header=None, names=[file_id])
                        else:
                            aux_data = pd.read_csv(file_path,sep='\t',index_col=0,header=None, names=[file_id])
                        # Concatenate previous results
                        if main_data is not None:
                            main_data = pd.concat([main_data, aux_data],axis=1,copy=False)
                        else:
                            main_data = aux_data
                        if i%1000 == 0:
                            print('{} files processed'.format(i))
                        i+=1
        # Calculate transposed matrix
        #main_data = main_data.transpose()
        # Write output file
        main_data.to_csv(output_file, index=True, sep=',')
    else:
        main_data = pd.read_csv(output_file, sep=',', index_col=0)

    print("The final dataframe contains {} rows (genes) and {} columns (samples)".format(len(main_data), len(main_data.columns)))
    
    return



#######################
# SECONDARY FUNCTIONS #
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

if  __name__ == "__main__":
    main()

