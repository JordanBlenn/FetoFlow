bc_utils
========

Boundary condition helpers used to validate and normalise inlet/outlet
specifications for simulations.

Overview
--------

The ``bc_utils`` module provides utilities to validate and construct the
boundary condition dictionary used by the solver. The primary function is
:func:`generate_boundary_conditions` which accepts either numeric (single)
or mapping (per-node/element) inlet specifications and a single outlet
pressure value. It performs type and value checks and returns a
normalized dict consumed by the matrix construction code.

generate_boundary_conditions
----------------------------

Create a normalized boundary condition dictionary used by the matrix
builders and solvers.

Parameters
----------
inlet_pressure : float or int or dict, optional
    If a single numeric value is provided it is used as the inlet
    pressure (Pa) applied at all inlet nodes. If a ``dict`` is supplied
    it should map 1-based node indices to numeric pressures and each
    value is validated and converted to internal 0-based indexing.

inlet_flow : float or int or dict, optional
    If a single numeric value is provided it is used as the total inlet
    flow. If a ``dict`` is supplied it should map 1-based element
    indices to flow values. Only one of ``inlet_pressure`` or
    ``inlet_flow`` may be supplied.

outlet_pressure : float or int, optional
    Outlet pressure (Pa) applied to terminal nodes. If omitted and an
    inlet flow is supplied a warning is emitted and a default of 0 is
    used.

Returns
-------
dict
    Normalised boundary conditions in the format::

        { 'inlet': {'pressure' or 'flow': value_or_dict}, 'outlet': {'pressure': value} }

Notes
-----

The returned mapping is the canonical format expected by the
:mod:`FetoFlow.matrix_builder` and :mod:`FetoFlow.pressure_flow_utils`
modules; prefer calling this function rather than assembling the dict
manually.

Examples
--------

Single inlet pressure::

    bc = generate_boundary_conditions(inlet_pressure=12000, outlet_pressure=0)

Multiple inlet pressures (per-node)::

    bc = generate_boundary_conditions(inlet_pressure={1: 12000, 5: 11000}, outlet_pressure=0)


API reference
-------------

.. automodule:: FetoFlow.bc_utils
    :members:
    :undoc-members:
    :show-inheritance:

