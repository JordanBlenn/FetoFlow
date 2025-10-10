matrix_builder
==============

Matrix builders for assembling the sparse linear systems used to solve
for pressures and flows. Two construction styles are provided:

- A full node/edge matrix builder (:func:`create_matrices`) that assembles
    a mixed system coupling nodal mass-conservation equations with edge
    pressure-drop equations. This produces a (n+m) x (n+m) sparse system
    where n is the number of nodes and m the number of edges.
- A reduced (``small``) matrix builder (:func:`create_small_matrices`)
    that uses network incidence matrices and a weighted Laplacian-style
    construction to produce a smaller r x r system. This is used for
    faster solves in typical tree-like topologies.

Important functions
-------------------

create_matrices(G, n, m, bcs)
        Build the (n + m) sparse linear system A and vector b. Performs
        validation of boundary conditions and constructs rows for:
        - mass conservation at nodes
        - pressure drop along edges (Ohm's law: p_u - p_v = R q)

create_small_matrices(G, bcs, branching_angles=False, non_linear_rheology=False)
        Build a reduced sparse system using network incidence matrices and a
        diagonal resistance weighting. Returns (A, b, bc_export, iter_options)
        where ``bc_export`` summarises boundary indices and ``iter_options``
        contains auxiliary matrices for iterative updates (branching-angle
        corrections etc.).

Notes
-----

- Boundary conditions for inlet/outlet must be created using
    :func:`FetoFlow.bc_utils.generate_boundary_conditions`.
- The reduced builder supports both pressure and flow inlet types and
    sets up the linear system accordingly.

API reference
-------------

.. automodule:: FetoFlow.matrix_builder
        :members:
        :undoc-members:
        :show-inheritance:

Function arguments
------------------

create_matrices
    Parameters
    ----------
    G : networkx.DiGraph
        Graph representing the vascular network. Edges must contain a
        ``resistance`` attribute and ``edge_id`` values.

    n : int
        Number of nodes in the graph (expected to equal ``G.number_of_nodes()``).

    m : int
        Number of edges in the graph (expected to equal ``G.number_of_edges()``).

    bcs : dict
        Boundary condition dictionary created by
        :func:`FetoFlow.bc_utils.generate_boundary_conditions`.

    Returns
    -------
    A : scipy.sparse.csr_matrix
        The assembled sparse system matrix of shape (n+m, n+m).

    b : numpy.ndarray
        Right-hand side vector of length (n+m).

        Notes
        -----
        - The function validates that the provided ``n`` and ``m`` match the
            graph size and will raise :class:`AssertionError` if they differ.
        - Use :func:`FetoFlow.bc_utils.generate_boundary_conditions` to
            prepare the ``bcs`` argument.

        Examples
        --------
        >>> A, b = create_matrices(G, n=G.number_of_nodes(), m=G.number_of_edges(), bcs=bcs)

create_small_matrices
    Parameters
    ----------
    G : networkx.DiGraph
        Graph representing the vascular network.

    bcs : dict
        Normalised boundary condition dict.

    branching_angles : bool, optional
        Include extra matrices to account for branching-angle pressure
        losses when returning iteration options.

    non_linear_rheology : bool, optional
        Reserved for future use; currently not fully implemented.

    Returns
    -------
    A : scipy.sparse matrix
        Reduced system matrix.

    b : numpy.ndarray
        Right-hand side vector for the reduced system.

    bc_export : tuple
        Tuple summarising boundary condition structure for the small system.

    iter_options : dict or None
        Auxiliary matrices and metadata needed for iterative updates
        (branching-angle matrices etc.) or None.

Examples
--------

.. code-block:: python

    # Full matrix
    A, b = create_matrices(G, n=G.number_of_nodes(), m=G.number_of_edges(), bcs=bcs)

    # Reduced small matrices
    A_s, b_s, bc_export, iter_opts = create_small_matrices(G, bcs)

Cross references
----------------

- The reduced matrices are consumed by :func:`FetoFlow.solve_utils.solve_small_system`.

    Examples
    --------
    >>> A, b, bc_export, iter_options = create_small_matrices(G, bcs, branching_angles=True)
    >>> # iterative solvers use iter_options to update A during non-linear solves

    See also
    --------
    :func:`FetoFlow.solve_utils.solve_small_system`, :func:`FetoFlow.solve_utils.iterative_solve_small`
