geometry_utils
==============

Graph-based geometry helpers. The package uses :mod:`networkx` directed
graphs to represent the vascular tree. Each edge stores attributes such as
``length``, ``radius``, ``resistance``, ``edge_id`` and ``strahler``.

Main responsibilities
---------------------

- Building the initial directed graph from node and element lists
    (:func:`create_geometry`).
- Adding a venous mesh based on arterial topology
    (:func:`create_venous_mesh`).
- Computing Strahler ordering for radii assignment
    (:func:`update_strahlers`).
- Utility helpers to compute lengths, branching angles and to add
    anastomoses (:func:`calcLength`, :func:`calculate_branching_angles`,
    :func:`create_anastomosis`).

create_geometry(nodes, elements, inlet_radius, strahler_ratio_arteries, ...)
        Construct a :class:`networkx.DiGraph` from raw node and element data,
        set per-edge defaults, compute Strahler ordering and assign radii.

        Parameters
        ----------
        nodes : dict
                Mapping ``node_id -> [x, y, z]`` (0-based indices). Produced by
                :func:`file_parsing_utils.read_nodes`.

        elements : list of tuple
                Connectivity list ``[(node_from, node_to), ...]`` (0-based indices).

        inlet_radius : float
                Radius (m) used for the root/inlet arterial vessels.

        strahler_ratio_arteries : float
                Scaling factor between Strahler orders used to set vessel radii.

        arteries_only : bool, optional
                When True do not create a venous mesh; only arterial graph is
                returned. Default is False.

        outlet_vein_radius : float, optional
                Outlet vein radius (m) used when building the venous mesh. Required
                if ``arteries_only`` is False.

        strahler_ratio_veins : float, optional
                Strahler scaling factor used for veins. Required when building
                venous mesh.

        fields : dict, optional
                Optional per-edge fields (for example radius or resistance) as
                produced by :func:`file_parsing_utils.define_fields_from_files`.

        default_mu : float, optional
                Default dynamic viscosity (Pa.s) assigned to edges.

        default_hematocrit : float, optional
                Default hematocrit value assigned to edges.

        Returns
        -------
        networkx.DiGraph
                Directed graph with nodes and edges annotated with attributes
                such as ``length``, ``radius``, ``edge_id``, ``strahler`` and
                ``vessel_type``.

        See also
        --------
        :func:`create_venous_mesh`, :func:`update_strahlers`, :func:`create_anastomosis`

        Examples
        --------
        >>> from FetoFlow.file_parsing_utils import read_nodes, read_elements
        >>> nodes = read_nodes('placenta.ipnode')
        >>> elements = read_elements('placenta.ipelem')
        >>> G = create_geometry(nodes, elements, inlet_radius=0.002, strahler_ratio_arteries=0.79)

create_venous_mesh(...)
        Build a venous graph by copying and relabelling the arterial graph,
        applying vein radii and reversing edge directions.

update_strahlers(G, node_in, node_out)
        Recursive Strahler calculation used during geometry construction.

create_anastomosis(G, node_from, node_to, radius=None, mu=...)
        Add a connecting edge (anastomosis) between two nodes and compute an
        approximate resistance for it. Input indices are in the external
        (1-based) IPNODE indexing and converted to internal 0-based indices.

Examples
--------

.. code-block:: python

        G = create_geometry(nodes, elements, inlet_radius=0.002, strahler_ratio_arteries=0.79)
        G = create_anastomosis(G, node_from=12, node_to=34, radius=0.1)

API reference
-------------

.. automodule:: FetoFlow.geometry_utils
        :members:
        :undoc-members:
        :show-inheritance:

Function arguments
------------------

create_geometry
        Parameters
        ----------
        nodes : dict
                Mapping ``node_id -> [x, y, z]`` (0-based indices). Produced by
                :func:`file_parsing_utils.read_nodes`.

        elements : list of tuple
                Connectivity list ``[(node_from, node_to), ...]`` (0-based indices).

        inlet_radius : float
                Radius (m) used for the root/inlet arterial vessels.

        strahler_ratio_arteries : float
                Scaling factor between Strahler orders used to set vessel radii.

        arteries_only : bool, optional
                When True do not create a venous mesh; only arterial graph is
                returned. Default is False.

        outlet_vein_radius : float, optional
                Outlet vein radius (m) used when building the venous mesh. Required
                if ``arteries_only`` is False.

        strahler_ratio_veins : float, optional
                Strahler scaling factor used for veins. Required when building
                venous mesh.

        fields : dict, optional
                Optional per-edge fields (for example radius or resistance) as
                produced by :func:`file_parsing_utils.define_fields_from_files`.

        default_mu : float, optional
                Default dynamic viscosity (Pa.s) assigned to edges.

        default_hematocrit : float, optional
                Default hematocrit value assigned to edges.

        Returns
        -------
        networkx.DiGraph
                Directed graph with nodes and edges annotated with attributes
                such as ``length``, ``radius``, ``edge_id``, ``strahler`` and
                ``vessel_type``.

                Notes
                -----
                The function will compute Strahler numbers for each edge and assign
                radii by scaling from the inlet radius. If a ``fields`` dict is
                provided and contains a ``radius`` field that value is used where
                specified and Strahler scaling applied elsewhere.

                Cross references
                ----------------
                - Use :func:`FetoFlow.file_parsing_utils.read_nodes` and
                  :func:`FetoFlow.file_parsing_utils.read_elements` to produce the
                  ``nodes`` and ``elements`` inputs to this function.

        create_venous_mesh
                Parameters
                ----------
                G : networkx.DiGraph
                        The arterial graph to copy and relabel into a venous mesh.

                Returns
                -------
                networkx.DiGraph
                        Venous mesh graph with relabelled nodes and vein attributes.

                Notes
                -----
                The venous mesh is generated by copying the arterial topology,
                relabelling node ids (offset by ``num_artery_nodes``) and reversing
                edges. The function is usually called internally by
                :func:`create_geometry` when ``arteries_only=False``.

        update_strahlers
                Parameters
                ----------
                G : networkx.DiGraph
                        Graph to update in-place.

                node_in, node_out : int
                        Edge endpoints describing the arc for which Strahler ordering
                        should be computed. This function is recursive and intended to
                        be used during geometry building.

                Returns
                -------
                networkx.DiGraph
                        The graph with updated ``strahler`` attributes.

        create_anastomosis
                Parameters
                ----------
                G : networkx.DiGraph
                        Graph to modify.

                node_from : int
                        1-based ipnode index identifying the anastomosis start node.

                node_to : int
                        1-based ipnode index identifying the anastomosis end node.

                radius : float, optional
                        Radius of the anastomosis (mm in calling code; will be converted to m).

                mu : float, optional
                        Viscosity used to compute an approximate resistance for the new edge.

                Returns
                -------
                networkx.DiGraph
                        The modified graph with the new anastomosis edge and attributes set.

                Examples
                --------

                .. code-block:: python

                        # build graph first from parsed files
                        G = create_geometry(nodes, elements, inlet_radius=0.002, strahler_ratio_arteries=0.79)
                        # add a small anastomosis between ipnode indices 12 and 34
                        G = create_anastomosis(G, node_from=12, node_to=34, radius=0.05)

create_venous_mesh
        Parameters
        ----------
        G : networkx.DiGraph
                The arterial graph to copy and relabel into a venous mesh.

        num_artery_nodes : int
                Number of arterial nodes (used for relabelling offsets).

        num_artery_edges : int
                Number of arterial edges (used to offset edge ids for veins).

        num_terminal_arterial_nodes : int
                Number of terminal arterial nodes (used when connecting capillary equivalents).

        outlet_vein_radius : float
                Vein outlet radius (m) used for radius assignment.

        strahler_ratio_veins : float
                Strahler ratio used to scale vein radii.

        max_strahler : int
                Maximum Strahler order in the arterial tree.

        Returns
        -------
        networkx.DiGraph
                Venous mesh graph with relabelled nodes and vein attributes.

create_anastomosis
        Parameters
        ----------
        G : networkx.DiGraph
                Graph to modify.

        node_from : int
                1-based ipnode index identifying the anastomosis start node.

        node_to : int
                1-based ipnode index identifying the anastomosis end node.

        radius : float, optional
                Radius of the anastomosis (mm in calling code; will be converted to m).

        mu : float, optional
                Viscosity used to compute an approximate resistance for the new edge.

        Returns
        -------
        networkx.DiGraph
                The modified graph with the new anastomosis edge and attributes set.
