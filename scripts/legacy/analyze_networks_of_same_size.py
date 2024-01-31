import sys, os
import graph_tool.all as gt
import networkx as nx
import numpy as np
import optparse
import pandas as pd
import random
import scipy
import time
import NetworkMedicineToolbox.wrappers as wrappers
import NetworkMedicineToolbox.network_utilities as network_utilities
import NetworkMedicineToolbox.parse_ncbi as parse_ncbi
import NetworkMedicineToolbox.parse_mesh as parse_mesh
import coexpression_analysis_functions_networkx as canx
random.seed(1510)


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_networks_of_same_size.py -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wto_N50 -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_sample_size -s 10 -r 10 -p 0.05
    """
    options = parse_options()
    analyze_gene_coexpression_networks_of_same_size(options)
    return


def parse_options():
    """
    This function parses the command line arguments and returns an optparse object.
    """

    parser = optparse.OptionParser("create_gene_coexpression_networks_cluster.py  -m METRIC")

    # Directory arguments
    parser.add_option("-n", action="store", type="string", dest="networks_dir", help="Directory with the input networks", metavar="NETWORKS_DIR")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Directory for the output files and plots", metavar="OUTPUT_DIR")
    parser.add_option("-g", action="store", type="string", dest="gene_expression_file", help="File containing the gene expression. It will be used to get the total number of genes", metavar="OUTPUT_DIR")
    parser.add_option("-s", action="store", type="string", dest="size", help="Sample size to analyze", metavar="SIZE")
    parser.add_option("-r", action="store", type="string", dest="repetitions", help="Number of repetitions to analyze", metavar="REPETITIONS")
    parser.add_option("-p", action="store", type="string", dest="pval_adj_cutoff", default=None, help="Cut-off of the adjusted p-value to consider an edge significant", metavar="PVAL_ADJ_CUTOFF")
    parser.add_option("-t", action="store", type="string", dest="top_percent", default=None, help="Top percentage of edges to select. If None, all of them are selected.", metavar="TOP_PERCENT")

    (options, args) = parser.parse_args()

    if options.networks_dir is None or options.output_dir is None or options.size is None or options.repetitions is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#################
#################
# MAIN FUNCTION #
#################
#################

def analyze_gene_coexpression_networks_of_same_size(options):
    """
    Analyze gene co-expression networks
    """

    # Define directories and files
    networks_dir = options.networks_dir
    output_dir = options.output_dir
    home_data_dir = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data'
    databases_dir = '/home/j.aguirreplans/Databases'
    ppi_dir = '/home/j.aguirreplans/data/PPI'

    # Parse all genes from gene expression data
    #gene_expression_file = '/home/j.aguirreplans/Projects/Scipher/SampleSize/data/bootstrap/nonresponders/RNAseq_NonResponders_all.txt'
    gene_expression_file = options.gene_expression_file
    #gene_expression_df = pd.read_csv(gene_expression_file, sep='\t')
    #genes_dataset = gene_expression_df.columns[1:]
    gene_expression_df = pd.read_csv(gene_expression_file, sep=",")
    genes_dataset = gene_expression_df.columns[1:]

    # Define parameters of bootstrap
    size = int(options.size)
    rep = int(options.repetitions)
    pval_adj_cutoff = options.pval_adj_cutoff
    if pval_adj_cutoff:
        pval_adj_cutoff = float(options.pval_adj_cutoff)
    top_percent = options.top_percent
    top_scoring = None
    if top_percent:
        num_nodes = len(genes_dataset)
        total_number_of_edges = float(num_nodes*(num_nodes-1))/2.0
        top_scoring = round(total_number_of_edges*float(top_percent)/100.0)
        print('Percentage of edges selected: {}% of {}. Number of edges selected: {}'.format(top_percent, int(total_number_of_edges), top_scoring))

    # Parse global network
    global_network_file = os.path.join(networks_dir, 'wto_RNAseq_NonResponders_all.net')
    #global_network_file = os.path.join(networks_dir, 'wto_RNAseq_Responders_all.net')
    #global_sig_network = canx.parse_significant_coexpression_network_from_wto_reading_the_file(global_network_file, pval_adj_cutoff=pval_adj_cutoff)
    global_sig_network = canx.parse_significant_coexpression_network_from_wto_with_pandas(global_network_file, pval_adj_cutoff=pval_adj_cutoff, top_scoring=top_scoring)

    print('Global gene co-expression network: {} nodes and {} edges'.format(len(global_sig_network.nodes()), len(global_sig_network.edges())))
    global_sig_network_edges = set(global_sig_network.edges())
    global_sig_network_nodes = set(global_sig_network.nodes())

    # Parse NCBI gene to get a dictionary gene ID to gene symbol
    ncbi_gene_info_file = os.path.join(home_data_dir, 'Homo_sapiens.gene_info')
    geneid_to_genesymbol, genesymbol_to_geneid = parse_ncbi.get_geneid_symbol_mapping(ncbi_gene_info_file)

    # Parse MeSH to get the class of each disease (the root concepts in MeSH)
    mesh_file = os.path.join(databases_dir, 'MeSH/mtrees2021.bin')
    m = parse_mesh.MESH(mesh_file)
    g = m.get_ontology(lower_concepts=True)

    # Parse disease genes associated to diseases with 20 genes or more
    # Map the disease to its corresponding class
    disease_genes_file = os.path.join(home_data_dir, 'Guney2016_GenesDisease.tsv')
    disease2genes = {}
    disease2genesdataset = {}
    disease2types = {}
    for i in open(disease_genes_file).readlines():
        v = i.rstrip().split('\t')
        disease = v[1]
        genes = v[2:]
        if len(genes) > 19:
            if disease in m.concept_to_concept_ids:
                disease2genes[disease] = [geneid_to_genesymbol[i] for i in genes if str(i) in geneid_to_genesymbol]
                disease2genesdataset[disease] = [gene for gene in disease2genes[disease] if gene in genes_dataset]
                concept_ids = m.concept_to_concept_ids[disease]
                for concept_id in concept_ids:
                    root_concept_id = concept_id.split('.')[0]
                    root_concept = m.concept_id_to_concept[root_concept_id]
                    disease2types.setdefault(disease, set()).add(root_concept)

    # Parse subsample networks
    topology_df = pd.DataFrame(columns=['size', 'rep', 'nodes', 'edges', 'av_degree', 'av_path_length', 'av_clustering_coef'])
    composition_df = pd.DataFrame(columns=['size', 'rep', 'nodes', 'edges', 'lost_or_gained'])
    components_df = pd.DataFrame(columns=['size', 'rep', 'num_components', 'size_lcc'])
    disease_df = pd.DataFrame(columns=['size', 'rep', 'disease', 'disease_class', 'total_disease_genes', 'total_disease_genes_in_dataset', 'disease_genes_in_network', 'lcc_size'])
    #kendall_df = pd.DataFrame(columns=['size', 'kendall'])
    wto_list_of_same_sample_size = []
    for r in range(1, rep+1):
        subsample_network_file = os.path.join(networks_dir, 'wto_RNAseq_NonResponders_sample_{}_rep_{}.net'.format(size, r))
        #subsample_network_file = os.path.join(networks_dir, 'wto_RNAseq_Responders_sample_{}_rep_{}.net'.format(size, r))
        if fileExist(subsample_network_file):

            start = time.time() # Starting time

            # Parse subsample network
            #subsample_sig_network = canx.parse_significant_coexpression_network_from_wto_reading_the_file(subsample_network_file, pval_adj_cutoff=pval_adj_cutoff)
            subsample_sig_network = canx.parse_significant_coexpression_network_from_wto_with_pandas(subsample_network_file, pval_adj_cutoff=pval_adj_cutoff, top_scoring=top_scoring)

            # Get weights & ranks associated to weights
            #wto_weights = [subsample_network[u][v]['wto'] for u,v in global_network.edges()]
            #wto_ranks = scipy.stats.rankdata(wto_weights)
            #wto_list_of_same_sample_size.append(wto_ranks)
            
            # Get significant network
            print('THE NETWORK HAS {} EDGES AND {} NODES'.format(len(subsample_sig_network.edges()), len(subsample_sig_network.nodes())))

            # Get components and LCC
            subsample_sig_components = network_utilities.get_connected_components(G=subsample_sig_network, return_as_graph_list=True)
            av_path_length_components = []
            for component in subsample_sig_components:
                # Calculate the average path length of each component (it cannot be calculated in disconnected graphs)
                av_path_length_components.append(nx.average_shortest_path_length(component))
            subsample_sig_lcc = subsample_sig_network.subgraph(subsample_sig_components[0])
            df2 = pd.DataFrame([[size, r, len(subsample_sig_components), len(subsample_sig_lcc.nodes())]], columns=['size', 'rep', 'num_components', 'size_lcc'])
            components_df = components_df.append(df2)
            subsample_sig_components = []
            subsample_sig_lcc = []

            # Get topological measures
            av_degree = 2*float(len(subsample_sig_network.edges()))/float(len(subsample_sig_network.nodes()))
            av_path_length = np.mean(av_path_length_components)
            av_clustering_coef = nx.average_clustering(subsample_sig_network)
            #communities = nx.algorithms.community.centrality.girvan_newman(subsample_sig_network)
            #node_groups = []
            #for com in communities:
            #  node_groups.append(list(com))
            #df2 = pd.DataFrame([[size, r, len(subsample_sig_network.nodes()), len(subsample_sig_network.edges()), av_degree, av_path_length, av_clustering_coef, len(node_groups)]], columns=['size', 'rep', 'nodes', 'edges', 'av_degree', 'av_path_length', 'av_clustering_coef', 'num_communities'])
            df2 = pd.DataFrame([[size, r, len(subsample_sig_network.nodes()), len(subsample_sig_network.edges()), av_degree, av_path_length, av_clustering_coef]], columns=['size', 'rep', 'nodes', 'edges', 'av_degree', 'av_path_length', 'av_clustering_coef'])
            topology_df = topology_df.append(df2)
            av_path_length_components = []
            communities = []
            #node_groups = []

            # Get connected disease genes
            for disease in disease2genes.keys():
                disease_genes = disease2genes[disease]
                disease_genes_dataset = disease2genesdataset[disease]
                disease_subgraph = subsample_sig_network.subgraph(disease_genes)
                disease_components = network_utilities.get_connected_components(G=disease_subgraph, return_as_graph_list=False)
                if len(disease_components) > 0:
                    disease_lcc = disease_subgraph.subgraph(disease_components[0])
                    disease_lcc_nodes = disease_lcc.nodes()
                else:
                    disease_lcc = []
                    disease_lcc_nodes = []
                for disease_class in disease2types[disease]:
                    df2 = pd.DataFrame([[size, r, disease, disease_class, len(disease_genes), len(disease_genes_dataset), len(disease_subgraph.nodes()), len(disease_lcc_nodes)]], columns=['size', 'rep', 'disease', 'disease_class', 'total_disease_genes', 'total_disease_genes_in_dataset', 'disease_genes_in_network', 'lcc_size'])
                    disease_df = disease_df.append(df2)
                disease_subgraph = []
                disease_components = []
                disease_lcc = []

            # Calculate differences between global and subsample network
            subsample_sig_network_edges = set(subsample_sig_network.edges())
            subsample_sig_network_nodes = set(subsample_sig_network.nodes())
            lost_nodes, gained_nodes, lost_edges, gained_edges = canx.calculate_difference_between_networks(global_sig_network_nodes, global_sig_network_edges, subsample_sig_network_nodes, subsample_sig_network_edges)
            #print(size, r, len(lost_edges), len(gained_edges))
            df3 = pd.DataFrame([[size, r, len(lost_nodes), len(lost_edges), 'lost']], columns=['size', 'rep', 'nodes', 'edges', 'lost_or_gained'])
            composition_df = composition_df.append(df3)
            df4 = pd.DataFrame([[size, r, len(gained_nodes), len(gained_edges), 'gained']], columns=['size', 'rep', 'nodes', 'edges', 'lost_or_gained'])
            composition_df = composition_df.append(df4)
            subsample_sig_network_edges = []
            subsample_sig_network_nodes = []
            lost_nodes = []
            lost_edges = []
            gained_nodes = []
            gained_edges = []

            end = time.time() # End time
            print("Runtime for sample size {} and rep {} is {}".format(size, r, end-start))

        else:
            print('Missing file: {}'.format(subsample_network_file))
        #break
    #W, count = calculate_kendall(np.array(wto_list_of_same_sample_size), n_permutations=1000)
    #df5 = pd.DataFrame([[size, W]], columns=['size', 'kendall'])
    #kendall_df = kendall_df.append(df5)

    topology_df.to_csv(os.path.join(output_dir, 'topology_size_{}.txt'.format(size)))
    components_df.to_csv(os.path.join(output_dir, 'components_size_{}.txt'.format(size)))
    disease_df.to_csv(os.path.join(output_dir, 'disease_components_size_{}.txt'.format(size)))
    composition_df.to_csv(os.path.join(output_dir, 'composition_size_{}.txt'.format(size)))
    #kendall_df.to_csv(os.path.join(output_dir, 'kendall_size_{}.txt'.format(size)))

    return

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