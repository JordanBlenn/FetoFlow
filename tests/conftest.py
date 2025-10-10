import sys
import pathlib

# Ensure src/ is on sys.path so tests can import the FetoFlow package directly.
ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pytest
import networkx as nx


@pytest.fixture
def small_linear_graph():
    """Create a small 3-node linear arterial graph for reuse in tests.

    Nodes: 0 -> 1 -> 2
    Edges: (0,1) edge_id=0, (1,2) edge_id=1
    """
    G = nx.DiGraph()
    for nid, coords in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]):
        G.add_node(nid, x=coords[0], y=coords[1], z=coords[2])

    edges = [(0, 1), (1, 2)]
    for eid, (u, v) in enumerate(edges):
        G.add_edge(
            u,
            v,
            edge_id=eid,
            resistance=1.0,
            length=1.0,
            radius=1e-3,
            strahler=1,
            vessel_type="artery",
            bifurcation_angle=0,
            mu=0.33600e-02,
            hematocrit=0.45,
            viscosity_factor=1,
        )
    return G


@pytest.fixture
def simple_two_node_graph():
    """Create a minimal two-node single-edge graph used by matrix_builder tests."""
    G = nx.DiGraph()
    G.add_node(0, x=0.0, y=0.0, z=0.0)
    G.add_node(1, x=1.0, y=0.0, z=0.0)
    G.add_edge(0, 1, edge_id=0, resistance=2.0, length=1.0, radius=1e-3, strahler=1, vessel_type="artery", mu=0.33600e-02, hematocrit=0.45, viscosity_factor=1)
    return G
