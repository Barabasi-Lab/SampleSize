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
    parser.add_option("-t", action="store", type="string", dest="type_gene_expression", default="tpm_unstranded", help="Type of gene expression (in case of parsing STAR files). Options: unstranded, stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded", metavar="TYPE_GENE_EXPRESSION")

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

    # Check type of gene expression
    type_gene_expression = options.type_gene_expression
    if type_gene_expression not in ["unstranded", "stranded_first", "stranded_second", "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"]:
        print("Unknown type of gene expression: {}".format(type_gene_expression))
        sys.exit(10)

    i=1
    if not fileExist(output_file):
        # Parse patients data
        main_data = None
        for patient_folder in folder_list:
            patient_folder_path = os.path.join(data_dir, os.fsdecode(patient_folder))
            if(os.path.isdir(patient_folder_path)):
                for file_id in os.listdir(patient_folder_path):
                    file_path = os.path.join(patient_folder_path, os.fsdecode(file_id))
                    if(len(file_id.split(".")) > 1 and (file_id.split(".")[-1] == "gz" or file_id.split(".")[-1] == "counts" or file_id.split(".")[-2] == "augmented_star_gene_counts")):
                        # Check type of gene expression file
                        if(file_id.split(".")[-1] == "gz" and len(file_id.split(".")) > 2 and file_id.split(".")[-3] == "htseq"):
                            f = gzip.open(file_path)
                            aux_data = pd.read_csv(f,sep='\t',index_col=0,header=None, names=[file_id])
                        elif(file_id.split(".")[-1] == "counts" and file_id.split(".")[-2] == "htseq"):
                            aux_data = pd.read_csv(file_path,sep='\t',index_col=0,header=None, names=[file_id])
                        elif(file_id.split(".")[-2] == "augmented_star_gene_counts"):
                            # Read data, skipping first line (irrelevant) and choosing first column (gene_id) as index column
                            aux_data = pd.read_csv(file_path,sep='\t',index_col=0, skiprows=1, header=0)
                            # Remove statistical columns
                            aux_data = aux_data[~aux_data.index.isin(["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"])]
                            # Select column with the type of gene expression selected and rename column using file_id
                            aux_data = aux_data[[type_gene_expression]].rename({type_gene_expression: file_id}, axis='columns')
                        else:
                            print("Unknown file format: {}".format(file_id))
                            sys.exit(10)
                        # Concatenate previous results
                        if main_data is not None:
                            #main_data = pd.concat([main_data, aux_data],axis=1,copy=False)
                            main_data = pd.merge(main_data, aux_data, left_index=True, right_index=True)
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

