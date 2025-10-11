resistance_utils
================

Compute hydraulic resistances for vessels and the capillary equivalent
elements. The module contains functions implementing simple Poiseuille
resistance for non-capillary vessels and a more detailed analytical
equivalent for capillary trees (the "analytical2015" model used by the
project).

Main functions
--------------

calculate_resistance(G, viscosity_model='constant', mu=..., capillary_model='analytical2015', capillary_parameters=None)
    Iterate over graph edges and set the ``resistance`` attribute on each
    edge. For non-capillary edges the code uses the Poiseuille formula
    8 mu L / (pi r^4). For capillary edges an equivalent resistance is
    computed using :func:`calculate_capillary_equivalent_resistance`.

calculate_viscosity_factor_from_radius(radius, hematocrit=0.45)
    Helper implementing a Pries-inspired empirical viscosity factor as a
    function of vessel radius and hematocrit.

calculate_capillary_equivalent_resistance(...)
    Compute the equivalent resistance of a tree/convolute capillary
    network given artery/vein inlet radii and capillary model
    parameters. This routine contains the analytical model used in the
    project and documents the expected parameters in ``capillary_parameters``.

Notes
-----

- The capillary model exposes many tunable parameters (number of
  series/parallel branches, convolute radius, lengths etc.). Default
  values are provided and validated by the calling code in
  :mod:`pressure_flow_utils`.

API reference
-------------

.. automodule:: FetoFlow.resistance_utils
    :members:
    :undoc-members:
    :show-inheritance:

Function arguments
------------------

calculate_resistance
    Parameters
    ----------
    G : networkx.DiGraph
        Graph whose edges will have their ``resistance`` attribute set.

    viscosity_model : str, optional
        One of ``'constant'``, ``'pries_network'``, ``'pries_vessel'`` or
        ``'flow_dependent'``. Controls how viscosity is handled for
        capillary and small vessels.

    mu : float, optional
        Base dynamic viscosity (Pa.s) used for Poiseuille calculations.

    capillary_model : str, optional
        Capillary equivalent model identifier (e.g. ``'analytical2015'``).

    capillary_parameters : dict, optional
        Parameters used by the capillary model. Keys include numbers of
        series/parallel convolutes, segment lengths, convolute radius,
        hematocrit and others. If omitted defaults are used.

    Returns
    -------
    networkx.DiGraph
        The modified graph with updated ``resistance`` attributes.

calculate_viscosity_factor_from_radius
    Parameters
    ----------
    radius : float
        Vessel radius in metres.

    hematocrit : float, optional
        Local hematocrit used by the empirical formula. Default ``0.45``.

    Returns
    -------
    float
        Empirical viscosity multiplication factor (dimensionless).

calculate_capillary_equivalent_resistance
    Parameters
    ----------
    radius_in_artery, radius_in_vein : float
        Radii (m) at the artery/vein sides used as boundary radii for the
        capillary equivalent calculation.

    viscosity_model, mu, capillary_model, capillary_parameters : see above

    Returns
    -------
    float
        Equivalent resistance (Pa s / m^3) of the capillary tree attached
        between the given arterial and venous radii.

Examples
--------

.. code-block:: python

    # Apply resistance calculation to a graph G
    G = calculate_resistance(G, viscosity_model='constant', mu=0.336e-2)

Cross references
----------------

- The capillary parameter dictionary used by :func:`calculate_capillary_equivalent_resistance` is assembled and validated in :func:`FetoFlow.pressure_flow_utils.pressures_and_flows`.

    Examples
    --------
    >>> G = calculate_resistance(G, viscosity_model='constant')

    See also
    --------
    :func:`FetoFlow.geometry_utils.create_anastomosis`, :func:`FetoFlow.pressure_flow_utils.pressures_and_flows`
