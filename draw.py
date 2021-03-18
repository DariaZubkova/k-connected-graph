import itertools
import networkx as nx
from typing import List, Tuple
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def init_graph(graphs): #-> nx.Graph:
    Graphs = []
    G = nx.Graph()
    Graph = []
    for graph in graphs:
        number = graph[0]
        nodes = graph[1]
        edges = graph[2]
        for i in range(len(nodes)):
            G.add_node(nodes[i])
        for i in range(len(edges)):
            G.add_edge(edges[i][0], edges[i][1], weight=edges[i][2])
        Graph.append(number)
        Graph.append(G)
        Graphs.append(Graph)
        G = nx.Graph()
        Graph = []
    return Graphs


def parse_stp(fname: str): #-> Tuple[list, List[tuple]]:
    nodes, edges, numbers = [], [], []
    number = 0
    edge = []
    graph = []
    graphs = []
    #fname = "C:\\Users\dasha\source\repos\graph\graph_k.txt"
    with open(fname, 'r') as inf:
        for line in inf.readlines():
            if len(line) < 2:
                if number != 0:
                    graph.append(number)
                    graph.append(nodes)
                    graph.append(edges)
                    graphs.append(graph)
                    nodes, edges, graph = [], [], []
                    number = 0
                continue
            if line[:2] == 'N ':
                _, n = line.split()
                number = n
                numbers.append(n)
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
        if number != 0:
            graph.append(number)
            graph.append(nodes)
            graph.append(edges)
            graphs.append(graph)
    return graphs

def draw(G):
    #pos = nx.spring_layout(G)
    #nx.draw(G, pos, node_color='#A0CBE2', edge_color='#BB0000', width=2, edge_cmap=plt.cm.Blues, with_labels=True)
    #nx.draw(G, pos, node_color='#A0CBE2', width=1, edge_cmap=plt.cm.Blues, with_labels=True)

    for graph in Graphs:
        #fig, ax = plt.subplot()
        fig = plt.figure(figsize=(7, 7)) #figsize=(9, 11)
        str = "Number graph is " + graph[0]
        plt.title(str)
        G = graph[1]
        nx.draw_circular(G, node_color='#A0CBE2', node_size=1000, with_labels=True)
        str_save = "graph" + graph[0] + ".png"
        #fig.tight_layout()
        plt.savefig(str_save, bbox_inches='tight') #, bbox_inches='tight'
        #plt.show()

    #nx.draw(G, with_labels=True)


def complete_graph(N: int) -> nx.Graph:
    graph = nx.Graph()

    N_range = range(N)

    all_pairs = itertools.permutations(N_range, 2)

    graph.add_nodes_from(N_range)
    graph.add_edges_from(all_pairs)

    return graph

if __name__ == "__main__":
    file_name = 'out64.txt'
    #nodes, edges, numbers
    graphs = parse_stp(file_name)
    Graphs = init_graph(graphs)
    #G = complete_graph(40)
    draw(Graphs)
    #plt.savefig("graph.png")
    #plt.show()