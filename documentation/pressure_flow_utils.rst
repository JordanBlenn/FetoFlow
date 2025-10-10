pressure_flow_utils
===================

High-level API used to run a full pressure-and-flow simulation from
placental mesh files and user parameters. This module ties together the
parsing, geometry construction, resistance calculation, matrix assembly
and solver steps.

Primary entry point
-------------------

pressures_and_flows(node_filename, element_filename, boundary_conditions, inlet_radius, strahler_ratio_arteries, *, ...)
        Run a complete pipeline:

        1. Read node/element files using :mod:`file_parsing_utils`.
        2. Validate and normalise :mod:`bc_utils` boundary conditions.
        3. Build vascular geometry with :mod:`geometry_utils`.
        4. Compute resistances (:mod:`resistance_utils`).
        5. Create sparse matrices (:mod:`matrix_builder`).
        6. Solve the system using :mod:`solve_utils` and export CSV outputs.

Key parameters
--------------

- ``node_filename`` / ``element_filename`` — input mesh files (``.ipnode``
    / ``.ipelem``) located in ``input_directory``.
- ``boundary_conditions`` — dict with possible keys ``inlet_pressure``,
    ``inlet_flow`` and ``outlet_pressure``. Use
    :func:`FetoFlow.bc_utils.generate_boundary_conditions` for convenience.
- ``inlet_radius`` / ``strahler_ratio_arteries`` — geometry parameters
    used to set arterial radii by Strahler ordering.
- ``viscosity_model`` — string selecting rheology handling. Supported
    options include ``constant`` and variants used by Pries-style models.

Outputs
-------

The function writes CSV files into ``output_directory`` containing node
pressures and per-element flows. Optionally it can return the graph with
pressures and flows saved as attributes by setting ``return_graph=True``.

Example
-------

.. code-block:: python

        bcs = { 'inlet_pressure': 12000, 'outlet_pressure': 0 }
        pressures_and_flows('placenta.ipnode', 'placenta.ipelem', bcs, inlet_radius=0.002, strahler_ratio_arteries=0.79)

API reference
-------------

.. automodule:: FetoFlow.pressure_flow_utils
        :members:
        :undoc-members:
        :show-inheritance:

Function arguments
------------------

pressures_and_flows
    Parameters
    ----------
    node_filename : str
        Name of the ``.ipnode`` file in ``input_directory``.

    element_filename : str
        Name of the ``.ipelem`` file in ``input_directory``.

    boundary_conditions : dict
        Dictionary that may contain keys ``inlet_pressure``, ``inlet_flow``
        and ``outlet_pressure``. See :func:`FetoFlow.bc_utils.generate_boundary_conditions`.

    inlet_radius : float
        Inlet radius (m) for arterial root vessels.

    strahler_ratio_arteries : float
        Strahler scaling factor used for arterial radius assignment.

    input_directory : str, optional
        Directory containing input files. Default is current directory (".").

    output_directory : str, optional
        Directory to write outputs. Default is ``./output_data``.

    flow_output_filename, pressure_output_filename : str, optional
        Filenames for CSV exports of flow and pressure.

    arteries_only : bool, optional
        If True do not generate a venous mesh; default False.

    viscosity_model : str, optional
        Viscosity handling mode. Supported values include ``constant``,
        ``pries_network``, ``pries_vessel`` and ``flow_dependent``.

    vessel_type : str, optional
        Placeholder for vessel elasticity model. Currently ``'rigid'`` is
        assumed.

    outlet_vein_radius, strahler_ratio_veins : float, optional
        Required when ``arteries_only`` is False; used to build the venous mesh.

    anastomosis : dict, optional
        Optional dict describing an anastomosis to add. Expected keys are
        ``node_from``, ``node_to`` and optional ``radius``.

    mu : float, optional
        Default dynamic viscosity used for non-capillary vessels (Pa.s).

    capillary_model : str, optional
        Capillary modelling choice (e.g. ``'analytical2015'``).

    capillary_parameters : dict, optional
        Parameters controlling the capillary equivalent model. See the
        source for defaults and expected keys.

    radius_filename : str, optional
        If provided, a field filename (``.ipfiel``) used to seed radii.

    other_field_filenames : dict, optional
        Mapping of additional field names to filenames to import.

    verbose, time_statistics, return_graph : bool, optional
        Misc flags controlling verbosity, timing output, and whether the
        function returns a graph object with pressures/flows attached.

    Returns
    -------
    networkx.DiGraph or None
        When ``return_graph=True`` the graph with ``pressure`` and
        ``flow`` attributes is returned. Otherwise the function writes
        CSV output and returns ``None``.

Examples and usage patterns
---------------------------

Single-file run producing CSV outputs

.. code-block:: python

    bcs = { 'inlet_pressure': 12000, 'outlet_pressure': 0 }
    pressures_and_flows('placenta.ipnode', 'placenta.ipelem', bcs, inlet_radius=0.002, strahler_ratio_arteries=0.79)

Return the graph for further inspection or plotting

.. code-block:: python

    G = pressures_and_flows('placenta.ipnode', 'placenta.ipelem', bcs, inlet_radius=0.002, strahler_ratio_arteries=0.79, return_graph=True)
    # Inspect per-edge flows
    flows = {d['edge_id']: d['flow'] for u,v,d in G.edges(data=True)}

Cross references
----------------

- Boundary dicts should be created using :func:`FetoFlow.bc_utils.generate_boundary_conditions`.
- Geometry is constructed with :func:`FetoFlow.geometry_utils.create_geometry` and resistances are calculated by :func:`FetoFlow.resistance_utils.calculate_resistance`.

        Notes
        -----
        - The function performs multiple validation steps on inputs and will
            raise :class:`ValueError` for incorrect types or missing required
            arguments (for example, ``outlet_vein_radius`` when building the
            venous mesh).
        - The default capillary parameters are used if ``capillary_parameters``
            is omitted; these defaults are warned about to the user.

        Examples
        --------
        Minimal run (arteries only)::

                bcs = {'inlet_pressure': 12000, 'outlet_pressure': 0}
                pressures_and_flows('placenta.ipnode', 'placenta.ipelem', bcs, inlet_radius=0.002, strahler_ratio_arteries=0.79, arteries_only=True)

        Full run with venous mesh and custom capillary parameters::

                cap_params = {'num_series':3, 'num_parallel':6, 'num_generations':3}
                bcs = {'inlet_pressure': 12000, 'outlet_pressure': 0}
                pressures_and_flows('n.ipnode', 'n.ipelem', bcs, inlet_radius=0.002, strahler_ratio_arteries=0.79, outlet_vein_radius=0.002, strahler_ratio_veins=0.8, capillary_parameters=cap_params)

        See also
        --------
        :func:`FetoFlow.geometry_utils.create_geometry`, :func:`FetoFlow.resistance_utils.calculate_resistance`, :func:`FetoFlow.solve_utils.solve_system`
