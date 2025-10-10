import pytest
import networkx as nx

from FetoFlow.pressure_flow_utils import pressures_and_flows


def _dummy_nodes():
    return {0: [0.0, 0.0, 0.0], 1: [1.0, 0.0, 0.0]}


def _dummy_elements():
    return [(0, 1)]


def _dummy_graph():
    G = nx.DiGraph()
    G.add_node(0, x=0.0, y=0.0, z=0.0)
    G.add_node(1, x=1.0, y=0.0, z=0.0)
    G.add_edge(0, 1, edge_id=0, resistance=1.0, length=1.0, radius=1e-3, strahler=1, vessel_type="artery", mu=0.33600e-02, hematocrit=0.45, viscosity_factor=1)
    return G


def test_invalid_file_extensions_raise():
    with pytest.raises(ValueError):
        pressures_and_flows(
            node_filename="nodes.txt",
            element_filename="elems.txt",
            boundary_conditions={},
            inlet_radius=1.0,
            strahler_ratio_arteries=0.5,
            input_directory=".",
            output_directory="./output_data",
        )


def test_missing_outlet_vein_radius_raises(monkeypatch, tmp_path):
    # monkeypatch file parsers to avoid needing real files
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_nodes", lambda p: _dummy_nodes())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_elements", lambda p: _dummy_elements())
    # create_geometry must succeed - return a small graph
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.create_geometry", lambda **kwargs: _dummy_graph())

    with pytest.raises(ValueError):
        pressures_and_flows(
            node_filename="nodes.ipnode",
            element_filename="elems.ipelem",
            boundary_conditions={"inlet_pressure": 100.0},
            inlet_radius=1.0,
            strahler_ratio_arteries=0.5,
            input_directory=str(tmp_path),
            output_directory=str(tmp_path / "out"),
            arteries_only=False,
        )


def test_invalid_inlet_radius_type_raises(monkeypatch, tmp_path):
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_nodes", lambda p: _dummy_nodes())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_elements", lambda p: _dummy_elements())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.create_geometry", lambda **kwargs: _dummy_graph())

    with pytest.raises(ValueError):
        pressures_and_flows(
            node_filename="nodes.ipnode",
            element_filename="elems.ipelem",
            boundary_conditions={"inlet_pressure": 100.0},
            inlet_radius="bad",
            strahler_ratio_arteries=0.5,
            input_directory=str(tmp_path),
            output_directory=str(tmp_path / "out"),
            arteries_only=True,
        )


def test_anastomosis_type_error(monkeypatch, tmp_path):
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_nodes", lambda p: _dummy_nodes())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_elements", lambda p: _dummy_elements())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.create_geometry", lambda **kwargs: _dummy_graph())

    with pytest.raises(TypeError):
        pressures_and_flows(
            node_filename="nodes.ipnode",
            element_filename="elems.ipelem",
            boundary_conditions={"inlet_pressure": 100.0},
            inlet_radius=1.0,
            strahler_ratio_arteries=0.5,
            input_directory=str(tmp_path),
            output_directory=str(tmp_path / "out"),
            arteries_only=True,
            anastomosis="not_a_dict",
        )


def test_invalid_capillary_model_raises(monkeypatch, tmp_path):
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_nodes", lambda p: _dummy_nodes())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.read_elements", lambda p: _dummy_elements())
    monkeypatch.setattr("FetoFlow.pressure_flow_utils.create_geometry", lambda **kwargs: _dummy_graph())

    with pytest.raises(ValueError):
        pressures_and_flows(
            node_filename="nodes.ipnode",
            element_filename="elems.ipelem",
            boundary_conditions={"inlet_pressure": 100.0},
            inlet_radius=1.0,
            strahler_ratio_arteries=0.5,
            input_directory=str(tmp_path),
            output_directory=str(tmp_path / "out"),
            arteries_only=True,
            capillary_model="not_supported",
        )
