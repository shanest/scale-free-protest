import numpy as np
import networkx as nx

# TODO: document the module, class, and experiment methods at end

class ProtestAgent(object):

    def __init__(self, active=False, threshold=2):
        self.active = active
        self._threshold = threshold

    @property
    def threshold(self):
        return self._threshold

    @property
    def active(self):
        return self._active

    @active.setter
    def active(self, boolean):
        self._active = boolean

    def activate(self):
        self.active = True


# https://stackoverflow.com/questions/28920824/generate-a-scale-free-network-with-a-power-law-degree-distributions
def scale_free_graph(num_nodes, gamma):
    """Generates a scale free graph of a certain size with a certain scale parameter.

    Args:
        num_nodes: size of the graph
        gamma: scaling parameter

    Returns:
        a networkx.Graph, which obeys a power-law with exponent gamma
    """
    scales = nx.utils.powerlaw_sequence(num_nodes, gamma)
    graph = nx.expected_degree_graph(scales, selfloops=False)
    return graph

def populate_graph(graph, threshold):
    """Populates a given graph with ProtestAgents.

    Args:
        graph: the graph to populate
        threshold: the threshold for the ProtestAgents

    Returns:
        a new graph, with graph.node now containing ProtestAgents
    """
    num_nodes = len(graph.nodes())
    graph.node = {i: ProtestAgent(threshold=threshold) for i in xrange(num_nodes)}
    #nx.set_node_attributes(graph, 'agent', agents_for_graph)
    return graph

def number_active_neighbors(graph, node):
    """Gets the number of active neighbors of a node in a graph.

    Args:
        graph: the graph
        node: the integer index of the node in the graph

    Returns:
        the number of ProtestAgent neighbors which are active
    """
    return np.sum([graph.node[neighbor_idx].active for neighbor_idx in graph[node].keys()])

def activate_nodes(graph, nodes, record_to=None):
    """Activate a given group of nodes in a graph.
    The arguments are modified, nothing is returned.

    Args:
        graph: the main graph
        nodes: the nodes to activate
        record_to: (optional) a set, containing nodes. If specified,
                    the newly activated nodes are unioned to it.
    """
    for agent in nodes:
        graph.node[agent].activate()
    if record_to is not None:
        # |= is union + assignment
        record_to |= set(nodes)

def run_trial(num_nodes, scaling_parameter, threshold, repression_rate):
    """Runs a trial of an experiment.  This method implements the basic logic of
    the spread of protest through a network, based on the number of an agent's
    neighbors who are already protesting.

    Args:
        num_nodes: the number of nodes in the graph to be built
        scaling_parameter: scale parameter for the power law that the graph will obey
        threshold: how many neighbors need to be protesting for an agent to begin protesting
        repression_rate: rate of node removal at each time step

    Returns:
        initial size: number of initially activated nodes
        initial density: density of initially activated subgraph
        initial clustering: average clustering coefficient of initially activated subgraph
        final size: number of protesting nodes at stop time
        num iters: how many iterations it took before stop condition was reached
    """
    graph = scale_free_graph(num_nodes, scaling_parameter)
    graph = populate_graph(graph, threshold)

    # INITIALIZE
    active_nodes = set([])
    # TODO: require seed_node to be in largest connected component?
    seed_node = np.random.randint(num_nodes)
    nodes_to_activate = [seed_node]
    nodes_to_activate.extend(graph[seed_node].keys())
    activate_nodes(graph, nodes_to_activate, active_nodes)

    # record some info
    initial_neighborhood = graph.subgraph(nodes_to_activate)
    initial_size = len(initial_neighborhood.nodes())
    initial_density = nx.density(initial_neighborhood)
    initial_clustering = nx.average_clustering(initial_neighborhood)

    # get ready
    num_iters = 0
    stop = False

    # MAIN LOOP
    # TODO: modify to incorporate repression for experiment 3
    while not stop:
        nodes_to_activate = []
        print len(active_nodes)

        # TODO: optimize this loop; only nodes added in previous step?
        for node in active_nodes:
            for neighbor in graph[node].keys():
                if (number_active_neighbors(graph, neighbor) >= graph.node[neighbor].threshold
                        and not graph.node[neighbor].active):
                    nodes_to_activate.append(neighbor)

        if nodes_to_activate == []:
            stop = True
        else:
            num_iters += 1
            activate_nodes(graph, nodes_to_activate, active_nodes)

    return initial_size, initial_density, initial_clustering, len(active_nodes), num_iters

def run_experiment(out_file, scales, repression_rates,
        num_nodes=40000, threshold=2, trials_per_setting=2500):
    """Runs an experiment.  Handles the main loops for running individual trials, as well
    as the recording of data to a file. Returns nothing, but writes to out_file.

    Args:
        out_file: file to write to
        scales: an iterable of possible scaling parameters
        repression_rates: an iterable of possible repression rates
        num_nodes: how many nodes to put in each graph
        threshold: the threshold to use for ProtestAgents
        trials_per_setting: how many trials to run per (scale X repression_rate) setting
    """
    data = []
    for gamma in scales:
        for repression_rate in repression_rates:
            parameters = (num_nodes, gamma, threshold, repression_rate)
            for _ in xrange(trials_per_setting):
                data.append(parameters + run_trial(num_nodes, gamma, threshold, repression_rate))

    head_line = ('num_nodes,gamma,threshold,repression_rate,initial_size,initial_density,' +
            'initial_clustering,final_size,num_iters')
    np.savetxt(out_file, data, delimiter=',', header=head_line, comments='')

def experiment_one(out_file='/tmp/exp1.csv'):
    """Runs experiment one, where no parameters vary.

    Args:
        out_file: file to write data to
    """
    #TODO: update to 2500
    run_experiment(out_file, [2.3], [0], trials_per_setting=2)

def experiment_two(out_file='/tmp/exp2.csv'):
    """Runs experiment two, where scale parameter varies.

    Args:
        out_file: file to write data to
    """
    scale_params = np.linspace(1, 3, num=200)
    run_experiment(out_file, scale_params, [0])

def experiment_three(out_file='/tmp/exp3.csv'):
    """Runs experiment three, where scale parameter and repression rate vary.

    Args:
        out_file: file to write data to
    """
    scale_params = np.linspace(1, 3, num=200)
    repression_rates = np.linspace(.1, 3, num=290)
    run_experiment(out_file, scale_params, repression_rates)

experiment_one()
