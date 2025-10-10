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

