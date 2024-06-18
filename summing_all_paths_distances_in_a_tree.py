from collections import defaultdict


# O(|V| + |E|) time complexity
def sumOfDistancesInTree(n: int, edges: list[list[int]]) -> list[int]:
    if n == 0:
        return []
    if n < 3:
        return [n - 1] * n

    graph = defaultdict(list)
    children_below = [0] * n
    children_dist_sum_below = [0] * n
    children_dist_sum_of = [0] * n

    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])

    def dfs_bottom_up(index: int, predecessor: int) -> None:
        neighbors = graph[index]

        if len(neighbors) == 1 and neighbors[0] == predecessor:
            children_below[index] = 0
            children_dist_sum_below[index] = 0
            return

        for neighbor in graph[index]:
            if neighbor == predecessor:
                continue
            dfs_bottom_up(neighbor, index)
            children_below[index] += children_below[neighbor] + 1
            children_dist_sum_below[index] += children_dist_sum_below[neighbor] + children_below[neighbor] + 1

    dfs_bottom_up(0, 0)

    def dfs_top_down(index: int, predecessor: int) -> None:
        neighbors = graph[index]

        if index == predecessor:
            children_dist_sum_of[index] = children_dist_sum_below[index]

        for neighbor in neighbors:
            if neighbor == predecessor:
                continue
            # n - 1 equals number of children of each node
            children_dist_sum_of[neighbor] = children_dist_sum_of[index] - children_below[neighbor] + \
                                             ((n - 1) - (children_below[neighbor] + 1))
            dfs_top_down(neighbor, index)

    dfs_top_down(0, 0)
    return children_dist_sum_of
