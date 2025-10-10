import numpy as np
import pytest

from FetoFlow.solve_utils import solve_small_system, solve_system


def test_solve_system_basic():
    # trivial solve for 2x2 linear system
    A = np.array([[2.0, 0.0], [0.0, 3.0]])
    b = np.array([4.0, 9.0])
    pressures, flows = solve_system(A=A, b=b, num_nodes=2, num_edges=0)
    assert pressures[0] == 2.0 and pressures[1] == 3.0


def test_solve_small_system_flow_and_pressure(simple_two_node_graph):
    G = simple_two_node_graph
    # Build a small Laplacian-like 1x1 system for pressure unknown at internal node
    # reuse create_small_matrices from matrix_builder in other tests; here build directly
    A = np.array([[1.0]])
    b = np.array([100.0])
    # Create a bc_export structure expected by solve_small_system: (type, indices, values, inlet_idx)
    bc_export = ("Pressure", np.array([0]), np.array([100.0]), None)
    pressures, flows = solve_small_system(A=A, b=b, G=G, boundary_conditions=bc_export)
    assert isinstance(pressures, dict)
    assert isinstance(flows, dict)

