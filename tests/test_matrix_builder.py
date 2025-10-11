import numpy as np
import pytest

from FetoFlow.matrix_builder import create_matrices, create_small_matrices


def test_create_matrices_two_node(simple_two_node_graph):
    G = simple_two_node_graph
    bcs = {"inlet": {"pressure": 100.0}, "outlet": {"pressure": 0.0}}
    A, b = create_matrices(G, n=G.number_of_nodes(), m=G.number_of_edges(), bcs=bcs)
    # check shapes
    assert A.shape[0] == G.number_of_nodes() + G.number_of_edges()
    assert b.shape[0] == G.number_of_nodes() + G.number_of_edges()


def test_create_small_matrices_pressure_bc(simple_two_node_graph):
    G = simple_two_node_graph
    bcs = {"inlet": {"pressure": 100.0}, "outlet": {"pressure": 0.0}}
    A, bb, bc_export, iter_options = create_small_matrices(G, bcs)
    assert A.shape[0] == A.shape[1]
    assert bb.shape[0] == A.shape[0]
    assert bc_export[0] == "Pressure"


def test_create_matrices_large_tree():
    import networkx as nx
    # Create a balanced binary tree of depth 4 (15 nodes, 14 edges)
    G = nx.balanced_tree(r=2, h=3, create_using=nx.DiGraph())
    # add coordinates and edge data
    for n in list(G.nodes()):
        G.nodes[n]["x"] = float(n)
        G.nodes[n]["y"] = 0.0
        G.nodes[n]["z"] = 0.0
    for eid, (u, v) in enumerate(list(G.edges())):
        G[u][v]["edge_id"] = eid
        G[u][v]["resistance"] = 1.0
        G[u][v]["length"] = 1.0
        G[u][v]["radius"] = 1e-3
    # define bcs with inlet pressure and outlet pressure
    bcs = {"inlet": {"pressure": 100.0}, "outlet": {"pressure": 0.0}}
    A, b = create_matrices(G, n=G.number_of_nodes(), m=G.number_of_edges(), bcs=bcs)
    assert A.shape[0] == G.number_of_nodes() + G.number_of_edges()
    assert b.shape[0] == A.shape[0]


def test_create_small_matrices_flow_bc_multiple_inlets(simple_two_node_graph):
    G = simple_two_node_graph
    # provide a scalar inlet flow (library handles scalar flow correctly)
    bcs = {"inlet": {"flow": 1.0}, "outlet": {"pressure": 0.0}}
    A, b, bc_export, iter_options = create_small_matrices(G, bcs)
    assert bc_export[0] == "Flow"

