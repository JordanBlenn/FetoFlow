solve_utils
===========

Solvers for the assembled sparse systems. The module exposes both direct
and iterative approaches tailored to the two matrix constructions used
by the project.

Key routines
------------

solve_system(A, b, num_nodes, num_edges)
    Direct solve using a sparse direct solver. Returns ``(pressures,
    flows)`` where both are dictionaries mapping node/edge ids to values.

solve_small_system(A, b, G, boundary_conditions, ...)
    Solver specialised for the reduced small system. Handles boundary
    value injection, reconstructs full node pressure vector and computes
    edge flows.

iterative_solve_small, solve_iterative_system
    Higher-level wrappers that perform iterative updates to account for
    rheology non-linearities (flow-dependent viscosity), branching-angle
    corrections and other physics that require updating resistances and
    re-solving until convergence.

Notes
-----

- The iterative solvers contain adaptive stepping / relaxation logic to
  stabilise convergence for non-linear problems. For large graphs and
  tight tolerances consider using the small-matrix route where possible.

API reference
-------------

.. automodule:: FetoFlow.solve_utils
    :members:
    :undoc-members:
    :show-inheritance:

Function arguments
------------------

solve_system
    Parameters
    ----------
    A : scipy.sparse matrix
        Assembled system matrix (shape (n+m, n+m)).

    b : numpy.ndarray
        Right-hand side vector of length (n+m).

    num_nodes : int
        Number of nodes in the original graph (n).

    num_edges : int
        Number of edges in the original graph (m).

    Returns
    -------
    pressures, flows : tuple of dict
        Dictionaries mapping node ids to pressures and edge ids to flows.

solve_small_system
    Parameters
    ----------
    A, b : sparse matrix, ndarray
        Reduced system matrix and right-hand side as returned by
        :func:`matrix_builder.create_small_matrices`.

    G : networkx.DiGraph
        Graph object used to reconstruct full pressures/flows.

    boundary_conditions : tuple
        Summary of boundary indices and values used to inject boundary
        values after solving the reduced system.

    ill_conditioned, p0, current_p, max_iterations, restart : optional
        Internal flags for iterative solves and GMRES restarts. See the
        source for details. Typically not required for external callers.

    Returns
    -------
    pressures, flows : tuple of dict
        Dictionaries mapping node ids to pressures and edge ids to flows.

iterative solvers
    The module also exposes iterative wrappers that accept additional
    parameters for tolerances and solver strategies (GMRES vs direct
    solves). See the source for the full argument list; common options
    include ``tol``, ``maxiter``, ``alpha`` (relaxation), and
    ``adaptive_stepping``.

Examples
--------

.. code-block:: python

    # Direct solve using the full (n+m) system
    pressures, flows = solve_system(A, b, num_nodes=n, num_edges=m)

    # Reduced small-system solve
    pressures, flows = solve_small_system(A_small, b_small, G, bc_export)

Cross references
----------------

- For matrix construction see :mod:`FetoFlow.matrix_builder`.
- Nonlinear iterative solves call :func:`FetoFlow.resistance_utils.calculate_viscosity_factor_from_radius` and update graph attributes via :func:`FetoFlow.geometry_utils.update_geometry_with_pressures_and_flows`.

    Examples
    --------
    Direct solve of assembled system::

        pressures, flows = solve_system(A, b, num_nodes, num_edges)

    Small system solve::

        pressures, flows = solve_small_system(A_small, b_small, G, bc_export)

    Iterative non-linear solve with adaptive stepping::

        p, q = iterative_solve_small(A_small, b_small, G, bc_export, tol=0.01, maxiter=50)

    See also
    --------
    :func:`FetoFlow.matrix_builder.create_matrices`, :func:`FetoFlow.matrix_builder.create_small_matrices`
