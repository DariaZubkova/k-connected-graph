import itertools
import networkx as nx
from typing import List, Tuple
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def init_graph(nodes, edges) -> nx.Graph:
    G = nx.Graph()
    for i in range(len(nodes)):
        G.add_node(nodes[i])
    for i in range(len(edges)):
        G.add_edge(edges[i][0], edges[i][1], weight=edges[i][2])
    return G


def parse_stp(fname: str) -> Tuple[list, List[tuple]]:
    nodes, edges = [], []
    edge = []
    #fname = "C:\\Users\dasha\source\repos\graph\graph_k.txt"
    with open(fname, 'r') as inf:
        for line in inf.readlines():
            if len(line) < 2:
                continue
            if line[:2] == 'E ':
                # E <u> <v> <cost>
                _, u, v, cost = line.split()
                edge.append(u)
                edge.append(v)
                edge.append(cost)
                edges.append(edge)
                edge = []
            elif line[:2] == 'T ':
                _, i = line.split()
                nodes.append(i)

    return nodes, edges

def draw(G):
    #pos = nx.spring_layout(G)
    #nx.draw(G, pos, node_color='#A0CBE2', edge_color='#BB0000', width=2, edge_cmap=plt.cm.Blues, with_labels=True)
    #nx.draw(G, pos, node_color='#A0CBE2', width=1, edge_cmap=plt.cm.Blues, with_labels=True)

    nx.draw_circular(G, node_color='#A0CBE2', node_size=1000, with_labels=True)

    #nx.draw(G, with_labels=True)


def complete_graph(N: int) -> nx.Graph:
    graph = nx.Graph()

    N_range = range(N)

    all_pairs = itertools.permutations(N_range, 2)

    graph.add_nodes_from(N_range)
    graph.add_edges_from(all_pairs)

    return graph

if __name__ == "__main__":
    file_name = 'g77.txt'
    nodes, edges = parse_stp(file_name)
    G = init_graph(nodes, edges)
    #G = complete_graph(40)
    draw(G)
    plt.savefig("graph.png")
    plt.show()