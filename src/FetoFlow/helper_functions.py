import networkx as nx


def getRadii(G: nx.digraph):
    """
    Returns all radii of all vessels in the Arterial tree as a dictionary:
    (nodeIn, nodeOut) : radius
    """
    radii = {}
    for u, v in G.edges():
        radii[(u, v)] = G[u][v]["radius"]
    return radii


def getRadius(G: nx.digraph, nodes: tuple = None):
    if nodes is not None:
        u, v = nodes
        return G[u][v]["radius"]
    else:
        raise ("Enter node values!")


def getNode(G, nodeID):
    return G.nodes()[nodeID]


def getEdgeData(G, edge):
    u, v = edge
    return G.get_edge_data(u, v)


def getVesselLength(G, initialVessel: tuple):
    """
    Check this
    """
    u, v = initialVessel
    length = G[u][v]["length"]
    out = G.out_edges(v)
    if len(out) == 1:
        length += getVesselLength(G, (v, out[0]))
    else:
        return length


def getNumVessels(G):
    """
    Returns number of distinct vessels in graph
    """
    successors = nx.dfs_sucessors(G, source=1)  # check if 0?
    count = 0
    if successors.get(1) == 1:
        count += 1
    for __, nodesOut in successors.items():
        unique = len(nodesOut)
        if unique > 1:
            count += unique
    return count


def branching_analytics():
    pass
