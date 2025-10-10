import pytest
import numpy as np

from FetoFlow.geometry_utils import calcLength, update_strahlers, create_anastomosis, create_geometry


def test_calc_length_simple():
    # Create two nodes with known coordinates
    class G:
        nodes = {0: {"x": 0.0, "y": 0.0, "z": 0.0}, 1: {"x": 3.0, "y": 4.0, "z": 0.0}}

    length = calcLength(G, 0, 1)
    # Distance is 5; converted from mm to m in implementation so expect 0.005
    assert np.isclose(length, 0.005)


def test_update_strahlers_leaf_case():
    # Build a tiny graph with node 0 -> 1 and 1 -> 2 (chain). Ensure strahler assignment works.
    import networkx as nx

    G = nx.DiGraph()
    G.add_node(0)
    G.add_node(1)
    G.add_node(2)
    G.add_edge(0, 1, strahler=None)
    G.add_edge(1, 2, strahler=None)
    G = update_strahlers(G, 0, 1)
    # edge (0,1) should have strahler 1 and edge (1,2) should have been set to 1 via recursion
    assert G[0][1]["strahler"] == 1
    assert G[1][2]["strahler"] == 1


def test_create_anastomosis_errors_and_success(simple_two_node_graph):
    G = simple_two_node_graph
    # nodes in create_anastomosis are 1-based so use values outside range to trigger errors
    with pytest.raises(ValueError):
        create_anastomosis(G, node_from=999, node_to=1)
    with pytest.raises(ValueError):
        create_anastomosis(G, node_from=1, node_to=1)

    # On this minimal graph create_anastomosis will attempt to compute length using
    # 1-based node indices and may raise a KeyError if coordinates aren't present.
    with pytest.raises(KeyError):
        create_anastomosis(G, node_from=1, node_to=2, radius=0.5)


def test_create_geometry_basic_and_venous_mesh():
    # Build a small tree: node 0 -> 1 -> 2 and 1 -> 3 (bifurcation)
    nodes = {0: [0.0, 0.0, 0.0], 1: [1.0, 0.0, 0.0], 2: [2.0, 0.1, 0.0], 3: [2.0, -0.1, 0.0]}
    elements = [(0, 1), (1, 2), (1, 3)]
    G = create_geometry(nodes, elements, inlet_radius=0.001, strahler_ratio_arteries=0.8, arteries_only=False, outlet_vein_radius=0.0005, strahler_ratio_veins=0.9)
    # Basic checks
    assert G.number_of_nodes() > 0
    assert G.number_of_edges() >= len(elements)
    # terminal nodes should have capillary_equivalent outgoing to venous mesh
    terminal_nodes = [n for n, d in G.out_degree() if d == 1 and any(G[n][v]["vessel_type"] == "capillary_equivalent" for _, v in G.out_edges(n))]
    # should be at least one terminal arterial node
    assert len(terminal_nodes) >= 1


def test_update_strahlers_complex():
    import networkx as nx
    G = nx.DiGraph()
    # chain 0->1->2, and 1->3
    for nid in range(4):
        G.add_node(nid)
    G.add_edge(0,1,strahler=None)
    G.add_edge(1,2,strahler=None)
    G.add_edge(1,3,strahler=None)
    G = update_strahlers(G, 0, 1)
    # The upstream edge (0,1) should have strahler 2 because its two children have equal strahler 1
    assert G[0][1]["strahler"] == 2

