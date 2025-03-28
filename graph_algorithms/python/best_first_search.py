from common_models import *
from typing import Callable
import heapq

def reconstruct_path(ancestors, start, goal):
    path_size = len(ancestors)
    path = [None] * path_size
    for i in reversed(range(len(ancestors))):
        path[i] = goal
        goal = ancestors[goal]
        if goal == start:
            path[i-1] = start
            break
    return path


# any heuristic - specifically inconsistent (possibly admissible) heuristic
# to determine what next node to visit is the most promising
# best first search uses f_score(n) = g_score(n) + heuristic(n) where g(n) is the cost to reach the node
def best_first_search(graph: Graph, start: Node, goal: Node, heuristic: Callable[Node, int]) -> list[Node]:
    ancestors : dict(Node, Node) = {}
    minheap = [(heuristic(start) + 0, start)]

    g_score = {node: float('inf') for node in graph.nodes}
    g_score[start] = 0
    f_score = {node: float('inf') for node in graph.nodes}
    f_score[start] = heuristic(start)

    while minheap:
        node = heapq.heappop(minheap)[1]
        if node == goal:
            return reconstruct_path(ancestors, start, goal)
        
        # print(f"currently at node {node}")

        for neighbor, edge_cost in graph.adjacency_list[node]:

            g_score_step = g_score[node] + edge_cost
            if g_score[neighbor] > g_score_step:
                ancestors[neighbor] = node
                g_score[neighbor] = g_score_step
                f_score[neighbor] = g_score_step + heuristic(neighbor)
                
                # print(f"neighbor {neighbor}, score estimate", g_score_step)
                # alternatively check if we only need to update the f_score in case it is already on the heap:
                heapq.heappush(minheap, (f_score[neighbor], neighbor))

    return "did not find a path"

# consistent heuristic - boils down to modified Dijkstra's algorithm
# with new distance d'(x,y) = d(x,y) + h(x) - h(y)


# example test for debugging
def main():
    test_graph = Graph()
    test_graph.add_node(1)
    test_graph.add_node(2)
    test_graph.add_node(3)
    test_graph.add_node(4)
    test_graph.add_node(5)
    test_graph.add_node(6)
    test_graph.adjacency_list = {
        1: [(2, 1), (3, 1), (4, 10)],
        2: [(1, 1), (5, 1), (6, 2)],
        3: [(1, 1), (5, 0)],
        4: [(1, 1), (6, 1)],
        5: [(2, 1), (3, 5), (4, 0)],
        6: [(2, 1)]
    }

    print(best_first_search(test_graph, 1, 6, lambda node: 1))

    test_graph = Graph()
    test_graph.add_node(1)
    test_graph.add_node(2)
    test_graph.add_node(3)
    test_graph.add_node(4)
    test_graph.add_node(5)
    test_graph.add_node(6)
    test_graph.adjacency_list = {
        1: [(2, 2), (5, 5)],
        2: [(3, 4)],
        3: [(4, 0)],
        4: [],
        5: [(6, 0)],
        6: [(4, 4), (2, 1), (3, 0)]
    }

    print(best_first_search(test_graph, 1, 4, lambda node: 1))


if __name__ == "__main__":
    main()

