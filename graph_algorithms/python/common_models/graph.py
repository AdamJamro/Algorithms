from typing import TypeAlias

class Node:
    def __init__(self, value = 0):
        self.value = value

    def __eq__(self, other):
        return self.value == other.value
    
    def __gt__(self, other):
        return self.value > other.value
    
    def __lt__(self, other):
        return self.value < other.value

    def __str__(self):
        return str(self.value)


class Edge:
    def __init__(self, node1, node2, weight=1):
        self.node1 = node1
        self.node2 = node2
        self.weight = weight

    def __str__(self):
        return f'({self.node1}, {self.node2})'


AdjacencyList: TypeAlias    = dict[Node, list[(Node, int)]]
AdjacencyMatrix : TypeAlias = list[list[(int, int)]]
EdgeList : TypeAlias        = list[Edge]


class Graph:
    def __init__(self):
        self.nodes : List[Node] = []
        self.edges : List[Edge] = []
        self.adjacency_list : AdjacencyList = {}
        self.adjacency_matrix : AdjacencyMatrix = []
        self.edge_list : EdgeList = []

    def add_node(self, node):
        self.nodes.append(node)
        self.adjacency_list[node] = []

    def add_edge(self, edge):
        self.edges.append(edge)
        self.adjacency_list[edge[0]].append(edge[1])
        self.adjacency_list[edge[1]].append(edge[0])

    def adjacency_matrix_from_adjacency_list(self):
        for node in self.nodes:
            row = []
            for neighbor in self.nodes:
                if neighbor in self.adjacency_list[node]:
                    row.append(1)
                else:
                    row.append(0)
            self.adjacency_matrix.append(row)


    def adjacency_matrix_from_edges(self):
        for node in self.nodes:
            row = []
            for neighbor in self.nodes:
                if (node, neighbor) in self.edges:
                    row.append(1)
                else:
                    row.append(0)
            self.edge_matrix.append(row)

    def edge_list(self):
        for edge in self.edges:
            self.edge_list.append(edge)

    def __str__(self):
        return str(self.adjacency_list)