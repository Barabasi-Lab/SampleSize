import graph_tool.all as gt
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
    network_df = pd.read_csv(network_file).drop_duplicates() # Read PPI and remove duplicates
    network_df = network_df.assign(equalcols=network_df.proteinA == network_df.proteinB)
    network_df = network_df[network_df['equalcols']==False] # Remove loops
    network = gt.Graph(directed=False)
    network_ids = network.add_edge_list(network_df[['proteinA','proteinB']].values, hashed=True)
    return network, network_ids


def parse_coexpression_network_from_wto(network_file):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network = gt.Graph(directed=False)
    wto = network.new_edge_property("double")
    pval_adj = network.new_edge_property("double")

    network_df = pd.read_csv(network_file, sep=",")
    network.add_edge_list(network_df.values, hashed=True, eprops=[wto, pval_adj])

    return network, wto, pval_adj


def filter_network_by_pval_adj(network, network_pval_adj_prop, pval_adj_cutoff=0.05):
    """
    Create a new network from the significant edges of another network.
    """
    #significant_pval_adj = network.new_ep("bool")
    # for edge in network.iter_edges():
    #     e = network.edge(*edge)
    #     if network_pval_adj_prop[e] < pval_adj_cutoff:
    #         significant_pval_adj[e] = 1
    #     else:
    #         significant_pval_adj[e] = 0
    # filtered_network = gt.GraphView(network, efilt=significant_pval_adj)
    filtered_network = gt.GraphView(network, efilt=lambda e: network_pval_adj_prop[e] < pval_adj_cutoff)

    #filtered_edges = {edge for edge in network.edges() if network_pval_adj_prop[edge] < pval_adj_cutoff}
    #filtered_network = gt.Graph(directed=False)
    #filtered_network.add_edge_list(filtered_edges)

    return filtered_network


def parse_significant_coexpression_network_from_wto(network_file, pval_adj_cutoff=0.05):
    """
    Parse co-expression network created using wTO.
    Columns of the wTO file: Node.1,Node.2,wTO,pval,pval.adj
    """
    network_df = pd.read_csv(network_file)
    network_df = network_df[network_df["pval.adj"] < pval_adj_cutoff]
    network = gt.Graph(directed=False)
    wto = network.new_edge_property("double")
    pval_adj = network.new_edge_property("double")
    #network_ids = network.add_edge_list(network_df[["Node.1", "Node.2"]].values, hashed=True)
    network_ids = network.add_edge_list(network_df.values, hashed=True, eprops=[wto, pval_adj])

    return network, network_ids, wto, pval_adj


def calculate_difference_between_networks(network1, network2):
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

