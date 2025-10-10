FetoFlow
========

Overview
--------

FetoFlow is a small library for constructing simplified placental vascular
networks from simple mesh files, computing vascular resistances (including
an analytical capillary equivalent model) and solving for nodal pressures
and element flows. The code is organised around a :mod:`networkx` directed
graph representation where each edge stores attributes such as ``length``,
``radius``, ``resistance`` and ``edge_id``.

The documentation is structured into module-level pages describing:

- parsing utilities (``file_parsing_utils``)
- geometry construction (``geometry_utils``)
- resistance computations (``resistance_utils``)
- matrix assembly (``matrix_builder``)
- solver routines (``solve_utils``)
- high-level pipeline (``pressure_flow_utils``)
- developer tools (``tools``)

Quick start
-----------

A minimal example of running a simulation::

    from FetoFlow.file_parsing_utils import read_nodes, read_elements
    from FetoFlow.pressure_flow_utils import pressures_and_flows

    bcs = {'inlet_pressure': 12000, 'outlet_pressure': 0}
    pressures_and_flows('placenta.ipnode', 'placenta.ipelem', bcs, inlet_radius=0.002, strahler_ratio_arteries=0.79)

See the module pages for detailed descriptions of functions and keyword
arguments.
