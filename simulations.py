from __future__ import division, print_function
from builtins import zip
from builtins import range
from builtins import object
from multiprocessing import Pool
import itertools
import argparse
from collections import defaultdict
from typing import Dict

import numpy as np
import networkx as nx
import community  # louvain
import pandas as pd
import tqdm

import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

# TODO: document the module


class ThresholdType(object):

    FIXED = "fixed"
    UNIFORM = "uniform"
    NORMAL = "normal"


class GraphType(object):

    SCALEFREE = "scale_free_graph"
    WATTS_STROGATZ = "watts_strogatz_graph"
    POWERLAW_CLUSTER = "powerlaw_cluster_graph"


class RepressionType(object):

    EDGE_REMOVAL = "edge_removal"
    NODE_REMOVAL = "node_removal"


class ProtestAgent(object):
    """A simple class, defining an agent who can be active/protesting or not.
    """

    def __init__(self, active=False, threshold=0.2):
        """Build a new ProtestAgent.

        Args:
            active: boolean, whether the agent is initially protesting or not
            threshold: int, how many neighbors need to be protesting before
                        this agent protests
        """
        self._active = active
        self._threshold = threshold

    @property
    def threshold(self):
        """int: the threshold -- number of neighbors that must be protesting --
        before this agent protests.
        """
        return self._threshold

    @property
    def active(self):
        """boolean: whether or not the agent is currently active/proesting. """
        return self._active

    @active.setter
    def active(self, boolean):
        self._active = boolean

    def activate(self):
        """Activate this agent.  Sets self.active to True. """
        self.active = True


def scale_free_graph(num_nodes, gamma):
    """Generates a scale free graph of a certain size with a certain
    scale parameter.

    Args:
        num_nodes: size of the graph
        gamma: scaling parameter

    Returns:
        a networkx.Graph, which obeys a power-law with exponent gamma
    """

    """Method 1: configuration model.  Fails on create_degree_sequence
    with low exponents.
    """
    """
    # TODO: networkx 2.0 removed create_degree_sequence!
    #scales = nx.utils.create_degree_sequence(num_nodes,
            #nx.utils.powerlaw_sequence, exponent=gamma, max_tries=500)
    scales = nx.utils.powerlaw_sequence(num_nodes, gamma)
    graph = nx.configuration_model(scales)
    loops = graph.selfloop_edges()
    graph = nx.Graph(graph)
    graph.remove_edges_from(loops)
    """
    graph = None
    # TODO: figure out ZeroDivisionError here, 1.4 seems OK, lower not...
    while graph is None:
        try:
            scales = nx.utils.powerlaw_sequence(num_nodes, gamma)
            graph = nx.expected_degree_graph(scales, selfloops=False)
        except:
            pass
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    return nx.Graph(graph.subgraph(components[0]))


def populate_graph(graph, threshold_fn):
    """Populates a given graph with ProtestAgents.

    Args:
        graph: the graph to populate
        threshold_fn: the fn to call for each ProtestAgent's threshold

    Returns:
        a new graph, with graph.node now containing ProtestAgents
    """
    for i in graph.nodes():
        new_threshold = threshold_fn()
        graph.nodes[i]["agent"] = ProtestAgent(threshold=new_threshold)
    return graph


def number_active_neighbors(graph, node):
    """Gets the number of active neighbors of a node in a graph.

    Args:
        graph: the graph
        node: the integer index of the node in the graph

    Returns:
        the number of ProtestAgent neighbors which are active
    """
    return np.sum(
        [
            graph.nodes[neighbor_idx]["agent"].active
            for neighbor_idx in graph[node].keys()
        ]
    )


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
        graph.nodes[agent]["agent"].activate()
    if record_to is not None:
        # |= is union + assignment
        record_to |= set(nodes)


def should_be_active(graph, node):
    """Checks whether a node should be active.

    Args:
        graph: a graph
        node: node in graph

    Returns:
        True iff the number of active neighbors of node in graph exceeds node's
        threshold
    """
    return (
        number_active_neighbors(graph, node) / len(list(graph[node].keys()))
        >= graph.nodes[node]["agent"].threshold
    )


def repress_edge_removal(graph, active_nodes, repression_rate):
    """Implements repression as edge removal.

    Args:
        graph: the main graph
        active_nodes: set of nodes currently active in graph
        repression_rate: probability of removing edge from active node
    """
    for node in active_nodes:
        neighbors = list(graph[node].keys())
        remove_which = np.random.binomial(1, repression_rate, size=(len(neighbors)))
        for idx in range(len(neighbors)):
            if remove_which[idx]:
                graph.remove_edge(node, neighbors[idx])


def repress_node_removal(graph, active_nodes, repression_rate, centralities):
    """Implements repression as targeted node removal.  Removes
    `repression_rate` percentage of the active nodes, in proportion to their
    degree centrality.

    Args:
        graph: the main graph
        active_nodes: set of nodes currently active in graph
        repression_rate: what % of active nodes to remove
        centralities: {node: degree} dict of degree centralities
    """
    to_remove = set()
    active = list(active_nodes)  # order needed for weights to match
    # TODO: make this more modular?
    # prob propto: norm(degrees) + (1-norm(thresholds))
    degrees = np.array([centralities[node] for node in active_nodes])
    degrees /= max(degrees)
    thresholds = np.array(
        [graph.nodes[node]["agent"].threshold for node in active_nodes]
    )
    thresholds /= max(thresholds)
    thresholds = [1 - threshold for threshold in thresholds]
    combined = degrees + thresholds
    probs = combined / sum(combined)
    num_to_remove = int(repression_rate * len(active))
    to_remove = set(np.random.choice(active, num_to_remove, replace=False, p=probs))
    # only remove nodes at end so that probabilities are from the same time
    graph.remove_nodes_from(to_remove)
    active_nodes -= to_remove


# TODO: just eliminate this method?
def repress_node_removal_old(graph, active_nodes):
    """Implements repression as targeted node removal.  The probability that an
    active node gets removed is proportional to its share of the activated
    eges.

    Args:
        graph: the main graph
        active_nodes: set of nodes currently active in graph
    """
    # list_active = list(active_nodes)
    num_neighbors = {node: len(list(graph.neighbors(node))) for node in active_nodes}
    total_neighbors = sum(num_neighbors.values())
    to_remove = set()
    for node in active_nodes:
        if np.random.random() < num_neighbors[node] / total_neighbors:
            to_remove.add(node)
    # only remove nodes at end so that probabilities are from the same time
    graph.remove_nodes_from(to_remove)
    active_nodes -= to_remove


def communities_with_protesters(partition, active_nodes):
    """Returns number of communities in a graph with protesters.

    Args:
        partition: node -> community ID dict
        active_nodes: set of nodes who are protesting

    Returns:
        integer, number of communities with protesters
    """
    return len(set([partition[node] for node in active_nodes]))


def protesting_communities(partition, active_nodes):
    """Gets the communities with protesters, and the size of each.

    Args:
        partition: node -> community ID dict
        active_nodes: set of nodes who are protesting

    Returns:
        community ID -> int dict, where int is the size of comm
    """
    communities = defaultdict(int)
    for node in active_nodes:
        communities[partition[node]] += 1
    return communities


def mean_community_protest_size(comm_sizes, protesting):
    if len(protesting) == 0:
        return 0
    return sum(protesting[comm] for comm in protesting) / sum(
        comm_sizes[comm] for comm in protesting
    )


def run_trial(
    num_nodes=1000,
    graph_type=GraphType.SCALEFREE,
    repression_type=RepressionType.NODE_REMOVAL,
    threshold_type=ThresholdType.NORMAL,
    **kwargs,
):
    """Runs a trial of an experiment.  This method implements the basic logic of
    the spread of protest through a network, based on the number of an agent's
    neighbors who are already protesting.

    Args:
        num_nodes: the number of nodes in the graph to be built
        scaling_parameter: scale parameter for the power law
                            that the graph will obey
        threshold: how many neighbors need to be protesting for an agent
                    to begin protesting
        repression_rate: rate of node removal at each time step

    Returns:
        initial size: number of initially activated nodes
        initial density: density of initially activated subgraph
        initial clustering: average clustering coefficient of
                            initially activated subgraph
        total nodes: size of network's largest component
        final size: number of protesting nodes at stop time
        num_iters: how many iterations it took before stopping
    """
    # BUILD GRAPH
    if graph_type == GraphType.SCALEFREE:
        graph = scale_free_graph(num_nodes, kwargs["scaling_parameter"])
    elif graph_type == GraphType.WATTS_STROGATZ:
        graph = nx.watts_strogatz_graph(num_nodes, kwargs["k"], kwargs["p"])
    elif graph_type == GraphType.POWERLAW_CLUSTER:
        graph = nx.powerlaw_cluster_graph(num_nodes, kwargs["m"], kwargs["p"])

    # POPULATE GRAPH WITH AGENTS
    if threshold_type == ThresholdType.FIXED:

        def threshold_fn():
            return kwargs["threshold"]

    elif threshold_type == ThresholdType.UNIFORM:
        threshold_fn = np.random.random
    elif threshold_type == ThresholdType.NORMAL:

        def threshold_fn():
            # TODO: make these kwargs
            return max(0, np.random.normal(0.25, 0.122))

    graph = populate_graph(graph, threshold_fn)

    # some global properties
    total_nodes = len(graph.nodes())
    # dict: {node: degree}
    centralities = nx.degree_centrality(graph)
    eigen_centralities = nx.eigenvector_centrality_numpy(graph)
    # betweenness = nx.betweenness_centrality(graph)
    # dict: {node: community ID}
    partition = community.best_partition(graph)
    community_sizes = defaultdict(int)
    for node in partition:
        community_sizes[partition[node]] += 1

    # INITIALIZE
    active_nodes = set([])
    seed_node = np.random.choice(graph.nodes())
    nodes_to_activate = [seed_node]
    nodes_to_activate.extend(list(graph[seed_node].keys()))
    activate_nodes(graph, nodes_to_activate, active_nodes)

    # record some info
    initial_neighborhood = graph.subgraph(nodes_to_activate)
    initial_size = len(initial_neighborhood.nodes())
    initial_density = nx.density(initial_neighborhood)
    initial_neighborhood_clustering = nx.average_clustering(initial_neighborhood)
    # clustering of initial nodes inside main graph
    initial_nodes_clustering = nx.average_clustering(graph, nodes=nodes_to_activate)
    # seed_degree = graph.degree[seed_node]
    initial_degrees = [graph.degree[node] for node in nodes_to_activate]
    initial_mean_degree = sum(initial_degrees) / len(initial_degrees)
    initial_median_degree = np.median(initial_degrees)
    seed_eigen_centrality = eigen_centralities[seed_node]
    initial_mean_eigen = sum(
        [eigen_centralities[node] for node in nodes_to_activate]
    ) / len(nodes_to_activate)
    # threshold statistics
    initial_mean_threshold = sum(
        [graph.nodes[node]["agent"].threshold for node in initial_neighborhood]
    ) / len(nodes_to_activate)
    initial_neighbors = set(
        [neighbor for neighbor in graph[node] for node in initial_neighborhood]
    )
    initial_neighbors_mean_threshold = sum(
        [graph.nodes[node]["agent"].threshold for node in initial_neighbors]
    ) / len(initial_neighbors)
    # initial_mean_betweenness = sum([betweenness[node] for node in nodes_to_activate]) / len(nodes_to_activate)
    # initial global measures
    initial_global_clustering = nx.average_clustering(graph)
    # avg_shortest_path = nx.average_shortest_path_length(graph)

    # DEFINE REPRESSION
    if repression_type == RepressionType.NODE_REMOVAL:

        def repress(graph, active_nodes):
            repress_node_removal(
                graph, active_nodes, kwargs["repression_rate"], centralities
            )

    elif repression_type == RepressionType.EDGE_REMOVAL:

        def repress(graph, active_nodes):
            repress_edge_removal(graph, active_nodes, kwargs["repression_rate"])

    # initial repression
    repress(graph, active_nodes)

    # store initial information
    communities = protesting_communities(partition, active_nodes)
    initial_dict = {
        "initial_size": initial_size,
        "initial_density": initial_density,
        "initial_neighborhood_clustering": initial_neighborhood_clustering,
        "initial_nodes_clustering": initial_nodes_clustering,
        # 'seed_degree': seed_degree,
        "seed_eigen_centrality": seed_eigen_centrality,
        "initial_median_degree": initial_median_degree,
        "initial_mean_degree": initial_mean_degree,
        "initial_mean_eigen_centrality": initial_mean_eigen,
        "initial_mean_threshold": initial_mean_threshold,
        "initial_neighbors_mean_threshold": initial_neighbors_mean_threshold,
        # 'initial_mean_betweenness_centrality': initial_mean_betweenness,
        "num_communities": len(set(partition.values())),
        "total_nodes": total_nodes,
        "initial_global_clustering": initial_global_clustering,
        "active_nodes": len(active_nodes),
        "communities_with_protesters": len(communities),
        "mean_community_protest_percent": mean_community_protest_size(
            community_sizes, communities
        ),
        "time_step": 0,
    }
    data = pd.DataFrame(initial_dict, index=[0])

    # get ready
    num_iters = 0
    stop = False

    # MAIN LOOP
    while not stop:

        # deactivate nodes that no longer should be active
        # because repression removes edges, formerly active can be de-activated
        # by no longer surpassing their threshold
        to_deactivate = set([])
        for node in active_nodes:
            if not should_be_active(graph, node):
                graph.nodes[node]["agent"].active = False
                to_deactivate.add(node)

        # initial neighborhood nodes can never be deactivated
        to_deactivate -= set(initial_neighborhood.nodes())
        active_nodes -= to_deactivate

        # Get set of neighbors that could be activated
        nodes_to_activate = []
        neighbors_set = []

        for node in active_nodes:
            neighbors_set.extend(list(graph[node].keys()))
        # only activate new candidates that are not already active
        neighbors_set = set(neighbors_set) - active_nodes

        for neighbor in neighbors_set:
            if (
                should_be_active(graph, neighbor)
                and not graph.nodes[neighbor]["agent"].active
            ):
                nodes_to_activate.append(neighbor)

        if nodes_to_activate == []:
            stop = True
        else:
            num_iters += 1
            activate_nodes(graph, nodes_to_activate, active_nodes)
            # repression
            repress(graph, active_nodes)
            # get info
            communities = protesting_communities(partition, active_nodes)
            iter_dict = dict(initial_dict)
            iter_dict["active_nodes"] = len(active_nodes)
            iter_dict["time_step"] = num_iters
            # TODO: record other measures per time step here!
            iter_dict["communities_with_protesters"] = communities_with_protesters(
                partition, active_nodes
            )
            iter_dict["communities_with_protesters"] = len(communities)
            iter_dict["mean_community_protest_percent"] = mean_community_protest_size(
                community_sizes, communities
            )
            data = data.append(iter_dict, ignore_index=True)

    data["num_iters"] = num_iters
    return data


# TODO: document!
def product_of_dict_lists(dicts):
    return [dict(list(zip(dicts, x))) for x in itertools.product(*list(dicts.values()))]


def run_trial_from_kw(keywords):
    results = run_trial(**keywords)
    for kw in keywords:
        results[kw] = keywords[kw]
    return results


def run_experiment(out_root, num_procs, trials_per_setting, **kwargs):
    # TODO: UPDATE DOCS!
    """Runs an experiment.  Handles the main loops for running individual
    trials, as well as the recording of data to a file. Returns nothing,
    but writes to out_file.

    Args:
        out_file: file to write to
        scales: an iterable of possible scaling parameters
        repression_rates: an iterable of possible repression rates
        num_nodes: how many nodes to put in each graph
        threshold: the threshold to use for ProtestAgents
        trials_per_setting: how many trials to run per
                            (scale X repression_rate) setting
        num_procs: how many processes to spawn to run trials
    """
    param_dicts = product_of_dict_lists(kwargs)
    out_file = (
        out_root
        + "-".join(
            [
                "{}={}".format(key, value[0])
                for key, value in kwargs.items()
                if len(kwargs[key]) == 1
            ]
        )
        + ".csv"
    )

    # build list of parameters
    parameters = []
    for param_idx in range(len(param_dicts)):
        params = param_dicts[param_idx]
        params["param_idx"] = param_idx
        for trial_idx in range(trials_per_setting):
            t_params = dict(params)  # copy so that can vary trial number
            t_params["trial_idx"] = trial_idx
            parameters.append(t_params)
        # parameters = [params for _ in range(trials_per_setting)]
    # send work to pool, wrapped in a progress bar
    procs = Pool(num_procs)
    data = pd.concat(
        list(
            tqdm.tqdm(procs.imap(run_trial_from_kw, parameters), total=len(parameters))
        ),
        ignore_index=True,
    )
    # write output
    data.to_csv(out_file)


def parse_config(exp_config: str) -> Dict:
    with open(f"exps/{exp_config}.yml", "r") as exp_file:
        config = yaml.load(exp_file, Loader=Loader)
    # parameters
    for param in ["repression_rate", "scaling_parameter", "p"]:
        if param in config:
            scale = getattr(np, config[param]['scale'])(**config[param]['args'])
            # add special value to some
            if 'insert' in config[param]:
                scale = np.insert(scale, *config[param]['insert'])
            config[param] = scale
    return config


if __name__ == "__main__":

    # TODO: make all exp options command line?
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", help="which experiment to run", type=str)
    parser.add_argument(
        "--num_procs", help="how many processes to use", type=int, default=1
    )
    parser.add_argument("--out_dir", help="path to output", type=str, default="/tmp")
    args = parser.parse_args()

    config = parse_config(args.exp)

    out_root = f"{args.out_dir}/{args.exp}-"
    run_experiment(out_root, num_procs=args.num_procs, **config)
