import networkx as nx
import os
import pandas as pd

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


def parse_ppi_network(network_file):
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


def parse_coexpression_network_from_wto_reading_the_file(network_file):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network = nx.Graph()
    with open(network_file) as network_fd:
        first_line = network_fd.readline()
        for line in network_fd:
            node1, node2, wto, pval, pval_adj = line.split(',')
            sorted_edge = sorted([node1, node2])
            network.add_edge(sorted_edge[0], sorted_edge[1], wto=float(wto), pval_adj=float(pval_adj))
    return network


def parse_coexpression_network_from_wto_with_pandas(network_file):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network_df = pd.read_csv(network_file, sep=",")
    network = nx.from_pandas_edgelist(network_df, 'Node.1', 'Node.2', ["wTO", "pval.adj"])
    return network


def filter_network_by_pval_adj(network, pval_adj_cutoff=0.05):
    """
    Create a new network from the significant edges of another network.
    """
    filtered_edges = [(u,v) for u,v,weights in network.edges(data=True) if weights['pval_adj'] < pval_adj_cutoff]  
    filtered_network = nx.Graph()
    filtered_network.add_edges_from(filtered_edges)
    return filtered_network


def parse_significant_coexpression_network_from_wto_reading_the_file(network_file, pval_adj_cutoff=0.05):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network = nx.Graph()
    with open(network_file) as network_fd:
        first_line = network_fd.readline()
        for line in network_fd:
            node1, node2, wto, pval, pval_adj = line.split(',')
            if float(pval_adj) < pval_adj_cutoff:
                sorted_edge = sorted([node1, node2])
                network.add_edge(sorted_edge[0], sorted_edge[1])
    return network


def parse_significant_coexpression_network_from_wto_with_pandas(network_file, pval_adj_cutoff=0.05, top_scoring=None):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network_df = pd.read_csv(network_file, sep=",")
    if pval_adj_cutoff:
        # Filter by adjusted p-value
        network_df = network_df[network_df["pval.adj"] < pval_adj_cutoff]
    if top_scoring:
        # Sort by column 3 (index 2), corresponding to score (wTO, WGCNA, pearson, spearman...) and keep only the top scoring edges
        network_df = network_df.sort_values(by=network_df[2], ascending=False, key=abs).head(top_scoring)
    # Keep the filtered network
    network = nx.from_pandas_edgelist(network_df, 'Node.1', 'Node.2')
    return network


def calculate_difference_between_networks(network1_nodes, network1_edges, network2_nodes, network2_edges):
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

