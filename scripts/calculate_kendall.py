import sys, os
import graph_tool.all as gt
import matplotlib.pyplot as plt
import multiprocessing as mp
import networkx as nx
import numpy as np
import optparse
import pandas as pd
import pylab
import random
import seaborn as sns
import scipy
import time
import NetworkMedicineToolbox.wrappers as wrappers
import NetworkMedicineToolbox.network_utilities as network_utilities
random.seed(1510)


def main():
    """
    module load python/3.8.1
    python /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/analyze_networks_of_same_size.py -n /scratch/j.aguirreplans/Scipher/SampleSize/networks_nonresponders_wto_N50 -o /home/j.aguirreplans/Projects/Scipher/SampleSize/data/out/analysis_sample_size -s 10 -r 10
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
    parser.add_option("-s", action="store", type="string", dest="size", help="Sample size to analyze", metavar="SIZE")
    parser.add_option("-r", action="store", type="string", dest="repetitions", help="Number of repetitions to analyze", metavar="REPETITIONS")

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
    ppi_dir = '/home/j.aguirreplans/data/PPI'

    # Define parameters of bootstrap
    size = int(options.size)
    rep = int(options.repetitions)
    pval_adj_cutoff = 0.05

    # Parse global network
    global_network_file = os.path.join(networks_dir, 'wto_RNAseq_NonResponders_all.net')
    global_network = parse_coexpression_network_from_wto_with_networkx(global_network_file)
    global_network_edges = set(global_network.edges())

    # Parse subsample networks
    kendall_df = pd.DataFrame(columns=['size', 'kendall'])
    wto_list_of_same_sample_size = []
    for r in range(1, rep+1):
        subsample_network_file = os.path.join(networks_dir, 'wto_RNAseq_NonResponders_sample_{}_rep_{}.net'.format(size, r))
        if fileExist(subsample_network_file):

            # Parse subsample network
            subsample_network = parse_coexpression_network_from_wto_with_networkx(subsample_network_file)

            # Get weights & ranks associated to weights
            wto_weights = [subsample_network[u][v]['wto'] for u,v in global_network_edges]
            wto_ranks = scipy.stats.rankdata(wto_weights)
            wto_list_of_same_sample_size.append(wto_ranks)

        else:
            print('Missing file: {}'.format(subsample_network_file))
        #break
    W, count = calculate_kendall(np.array(wto_list_of_same_sample_size), n_permutations=1000)
    df5 = pd.DataFrame([[size, W]], columns=['size', 'kendall'])
    kendall_df = kendall_df.append(df5)
    kendall_df.to_csv(os.path.join(output_dir, 'kendall_size_{}.txt'.format(size)))

    return

def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)

def parse_coexpression_network_from_wto_with_networkx(network_file):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    #network_df = pd.read_csv(network_file, sep=",")
    #edges = zip(network_df["Node.1"], network_df["Node.2"], network_df["wTO"])
    network = nx.Graph()
    with open(network_file) as network_fd:
        first_line = network_fd.readline()
        for line in network_fd:
            node1, node2, wto, pval, pval_adj = line.split(',')
            sorted_edge = sorted([node1, node2])
            network.add_edge(sorted_edge[0], sorted_edge[1], wto=float(wto), pval_adj=float(pval_adj))
    return network

def parse_significant_coexpression_network_from_wto_with_networkx(network_file, pval_adj_cutoff=0.05):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    #network_df = pd.read_csv(network_file, sep=",")
    #edges = zip(network_df["Node.1"], network_df["Node.2"], network_df["wTO"])
    network = nx.Graph()
    with open(network_file) as network_fd:
        first_line = network_fd.readline()
        for line in network_fd:
            node1, node2, wto, pval, pval_adj = line.split(',')
            if float(pval_adj) < pval_adj_cutoff:
                sorted_edge = sorted([node1, node2])
                network.add_edge(sorted_edge[0], sorted_edge[1])
    return network

def parse_coexpression_network_from_wto_with_graphtools(network_file):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network = gt.Graph(directed=False)
    network_wto_prop = network.new_edge_property("double")
    network.edge_properties['wto'] = network_wto_prop
    network_pval_adj_prop = network.new_edge_property("double")
    network.edge_properties['pval_adj'] = network_pval_adj_prop
    #network_df = pd.read_csv(network_file, sep=",")
    #print(network_df.values)
    #network.add_edge_list(network_df.values, hashed=True, eprops=)
    with open(network_file) as network_fd:
        first_line = network_fd.readline()
        for line in network_fd:
            node1, node2, wto, pval, pval_adj = line.split(',')
            network.add_edge(node1, node2)
            edge = network.add_edge(node1, node2)
            network_wto_prop[edge] = float(wto)
            network_pval_adj_prop[edge] = float(pval_adj)
    return network, network_wto_prop, network_pval_adj_prop

def parse_ppi_network_with_networkx(network_file):
    """
    Parse protein-protein interaction network.
    Columns of the file: proteinA,proteinB,databases
    """
    ppi = pd.read_csv(network_file).drop_duplicates() # Read PPI and remove duplicates
    ppi = ppi.assign(equalcols=ppi.proteinA == ppi.proteinB)
    ppi = ppi[ppi['equalcols']==False] # Remove loops
    edges = zip(ppi.proteinA, ppi.proteinB)
    ppi = nx.Graph()
    ppi.add_edges_from(edges)
    return ppi

def parse_ppi_network_with_graphtools(network_file):
    """
    Parse protein-protein interaction network.
    Columns of the file: proteinA,proteinB,databases
    """
    network_df = pd.read_csv(network_file).drop_duplicates() # Read PPI and remove duplicates
    network_df = network_df.assign(equalcols=network_df.proteinA == network_df.proteinB)
    network_df = network_df[network_df['equalcols']==False] # Remove loops
    network = gt.Graph(directed=False)
    network.add_edge_list(network_df[['proteinA','proteinB']].values, hashed=True)
    return network

def filter_network_by_pval_adj_with_networkx(network, pval_adj_cutoff=0.05):
    """
    Create a new network from the significant edges of another network.
    """
    filtered_edges = [(u,v) for u,v,weights in network.edges(data=True) if weights['pval_adj'] < pval_adj_cutoff]  
    filtered_network = nx.Graph()
    filtered_network.add_edges_from(filtered_edges)
    return filtered_network

def filter_network_by_pval_adj_with_graphtools(network, network_pval_adj_prop, pval_adj_cutoff=0.05):
    """
    Create a new network from the significant edges of another network.
    """
    filtered_edges = {edge for edge in network.edges() if network_pval_adj_prop[edge] < pval_adj_cutoff}
    filtered_network = gt.Graph(directed=False)
    filtered_network.add_edge_list(filtered_edges)
    return filtered_network

def calculate_difference_between_networks_with_networkx(network1_nodes, network1_edges, network2_nodes, network2_edges):
    """
    Calculate the nodes and edges lost and gained between network1 and network2
    """
    #network1 = nx.Graph()
    #network1.add_edges_from([('a', 'b'), ('c', 'd')])
    #network2 = nx.Graph()
    #network2.add_edges_from([('a', 'b'), ('e', 'f')])
    #network1_nodes = set(network1.nodes())
    #network2_nodes = set(network2.nodes())
    lost_nodes = network1_nodes - network2_nodes
    gained_nodes = network2_nodes - network1_nodes
    #network1_edges = set([frozenset(edge) for edge in network1.edges()])
    #network2_edges = set([frozenset(edge) for edge in network2.edges()])
    lost_edges = network1_edges - network2_edges
    gained_edges = network2_edges - network1_edges
    #print(len(network1_edges), len(network2_edges), len(lost_edges), len(gained_edges))
    return lost_nodes, gained_nodes, lost_edges, gained_edges

def calculate_difference_between_networks_with_graphtools(network1, network2):
    """
    Calculate the nodes and edges lost and gained between network1 and network2
    """
    network1_nodes = set(network1.get_vertices())
    network2_nodes = set(network2.get_vertices())
    lost_nodes = network1_nodes - network2_nodes
    gained_nodes = network2_nodes - network1_nodes
    network1_edges = set([frozenset(edge) for edge in network1.get_edges()])
    network2_edges = set([frozenset(edge) for edge in network2.get_edges()])
    lost_edges = network1_edges - network2_edges
    gained_edges = network2_edges - network1_edges
    return lost_nodes, gained_nodes, lost_edges, gained_edges

def calculate_kendall(values, n_permutations=1000):
    """
    Code based on: https://stackoverflow.com/questions/48893689/kendalls-coefficient-of-concordance-w-in-python
    """
    m = values.shape[0]
    n = values.shape[1]
    W = kendall_w(values)
    count = 0
    for trial in range(n_permutations):
        perm_trial = []
        for _ in range(m):
            perm_trial.append(list(np.random.permutation(range(1, 1+n))))
        count += 1 if kendall_w(np.array(perm_trial)) > W else 0
    print ('Calculated value of W:', W, ' exceeds permutation values in', count, 'out of 1000 cases')
    return W, count

def kendall_w(expt_ratings):
    """
    Code from: https://stackoverflow.com/questions/48893689/kendalls-coefficient-of-concordance-w-in-python
    """
    if expt_ratings.ndim!=2:
        raise 'ratings matrix must be 2-dimensional'
    m = expt_ratings.shape[0] #raters
    n = expt_ratings.shape[1] # items rated
    denom = m**2*(n**3-n)
    rating_sums = np.sum(expt_ratings, axis=0)
    S = n*np.var(rating_sums)
    return 12*S/denom

if  __name__ == "__main__":
    main()