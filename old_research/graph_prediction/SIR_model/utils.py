# Network Analysis of Vaccination Strategies
# Copyright (C) 2020 by The RAND Corporation
# See LICENSE.txt and README.txt for information on usage and licensing

import collections
import networkx as nx
import numpy as np
import pandas as pd
import scipy



def mean_confidence_interval(data, confidence=0.95):
    '''given an array of data, return the mean and confidence intervals'''
    n = data.shape[0]
    m, se = np.mean(data, axis=0), scipy.stats.sem(data, axis=0)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, np.maximum(m-h, 0*(m-h)), m+h


def neighbor_graph(G, root, k=2):
    '''return the induced k-degree or k-hop subgraph'''
    edges = nx.bfs_edges(G, root, depth_limit=k)
    nodes = [root] + [v for u, v in edges]
    return G.subgraph(nodes)


def graph_from_edge_list(edge_list):
    '''
    Build a graph from given an attributed edge list.
    It is not assumed that each edge is unique.
    '''
    G = nx.Graph()
    weight_dict = {}
    processed_edges = set()

    for x in edge_list:

        x0, x1 = np.sort([int(x[0]), int(x[1])])
        weight = x[2]/3600/24 #use units of days
        if not (x0, x1) in processed_edges:
            G.add_edge(x0, x1, weight=weight)
            processed_edges.add((x0, x1))
            weight_dict[(x0, x1)] = weight
        else:
            G[x0][x1]['weight'] += weight
            weight_dict[(x0, x1)] += weight

    return G


def print_network_stats(G):
    '''print out some useful statistics about the graph'''

    print('number of edges %i ' %G.number_of_edges())
    print('number of nodes %i ' %G.number_of_nodes())
    print('number of edges/nodes %.2f ' %(G.number_of_edges()/G.number_of_nodes()))
    print('density %.2e' %nx.density(G))

    largest_cc = max(nx.connected_components(G), key=len)
    print('fraction of nodes in largest component: %.3f' %(len(largest_cc)/G.number_of_nodes()))

    return


def heterogeneity(G):
    '''compute the <k>, <k>^2, kappa'''

    nodes = list(G.nodes())
    degrees = [G.degree[nodes[i]] for i in range(len(nodes))]

    result = {'E(k)': np.mean(degrees),
              'E(k^2)': np.var(degrees) + np.mean(degrees)**2,
              'kappa': (np.var(degrees) + np.mean(degrees)**2)/np.mean(degrees),
              'kappa/E(k)': (np.var(degrees) + np.mean(degrees)**2)/np.mean(degrees)**2
             }

    return result


def degree_getter(G, k):
    '''return all the degree-k nodes'''
    tmp = []
    for i, j in list(G.degree()):
        if j == k:
            tmp.append(i)

    return np.asarray(tmp)


def degree_getter_range(G, kmin, kmax):
    '''return all the nodes with kmin < degree < kmax'''
    tmp = []
    for i, j in list(G.degree()):
        if kmin <= j and j <= kmax:
            tmp.append(i)

    return np.asarray(tmp)


def degree_getter_improved(G, k):
    '''
    Return a random degree-k node, and if there is no such node,
    return a random node with the next highest degree.
    '''

    nodes = degree_getter(G, k)
    if len(nodes) != 0:
        return np.random.choice(nodes)
    else:
        degrees = np.asarray(list(dict(G.degree).values()))
        degrees = degrees[degrees < k]
        return degree_getter_improved(G, np.max(degrees))


def degree_dist_calculator(G, nodes=None):
    '''
    Helper function that computes the degree counts for a given network.
    '''
    degree_sequence = np.asarray(sorted([d for n, d in G.degree(nodes)], reverse=False))
    degree_counter = collections.Counter(degree_sequence)
    return np.asarray(list(degree_counter.keys())), np.asarray(list(degree_counter.values()))


def interpolation(tlist, flist):

    ## use one of the time-step sequences for the final evaluation
    ## use the one which terminates the earliest
    imin = np.argmin([tlist[i][-1] for i in range(len(tlist))])
    t = tlist[imin]

    flist_interpolated = []
    for i in range(len(flist)):
        f_interpolated = scipy.interpolate.interp1d(tlist[i], flist[i])
        flist_interpolated.append(f_interpolated(t))

    return t, np.asarray(flist_interpolated)
