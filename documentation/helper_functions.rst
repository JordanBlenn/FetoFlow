helper_functions
================

Small convenience helpers for interrogating the networkx graph. These are
thin wrappers used by higher-level code and tests.

Key functions
-------------

- :func:`getRadii` — return a mapping of edge tuples to radii.
- :func:`getRadius` — get the radius for a particular edge tuple.
- :func:`getNode` — access node attributes.
- :func:`getEdgeData` — return the edge data dictionary.
- :func:`getVesselLength` — compute the length of a vessel by following
  connected edges recursively.

Examples
--------

.. code-block:: python

    radii = getRadii(G)
    r = getRadius(G, (3, 4))

API reference
-------------

.. automodule:: FetoFlow.helper_functions
    :members:
    :undoc-members:
    :show-inheritance:

Function arguments
------------------

getRadii
  Parameters
  ----------
  G : networkx.DiGraph
    Graph containing vessel edges annotated with a ``radius``
    attribute.

  Returns
  -------
  dict
    Mapping ``(u, v) -> radius`` for each edge.

getRadius
  Parameters
  ----------
  G : networkx.DiGraph
    The graph.

  nodes : tuple, optional
    Edge tuple ``(u, v)``. If omitted the function raises.

  Returns
  -------
  float
    Radius of the requested edge.

getVesselLength
  Parameters
  ----------
  G : networkx.DiGraph
    The graph.

  initialVessel : tuple
    Edge tuple ``(u, v)`` from which to accumulate connected lengths.

  Returns
  -------
  float
    Cumulative length of the vessel (m).

Examples
--------

.. code-block:: python

  radii = getRadii(G)
  r = getRadius(G, (3, 4))
  length = getVesselLength(G, (0, 1))

Cross references
----------------

- These helpers are used extensively by the higher-level solvers in
  :mod:`FetoFlow.solve_utils` and by geometry utilities in
  :mod:`FetoFlow.geometry_utils`.

  Examples
  --------
  >>> getVesselLength(G, (0, 1))
