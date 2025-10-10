from warnings import warn
def generate_boundary_conditions(inlet_pressure=None, inlet_flow=None, outlet_pressure=None):
    # TODO: IF ONLY INLET FLOW, SET OUTLET PRESSURE TO 0 AND WE LOOK AT THE PRESSURE DROP. THEN NOTE THIS SOMEWHERE.
    if inlet_pressure and inlet_flow:
        raise TypeError(f"Cannot have both an inlet pressure and inlet flow defined. Inlet pressure: {inlet_pressure}Inlet Flow: {inlet_flow}")
    elif inlet_pressure is None and inlet_flow is None:
        raise TypeError("No inlet boundary condition defined. Must define one of inlet_pressure or inlet flow for a valid boundary condition.")
    bcs = {}

    if outlet_pressure is None and inlet_flow is None:
        raise TypeError(
            "No valid outlet pressure for boundary condition. "
            "Currently Reprosim does not support outlet flow boundary conditions and must have a defined outlet pressure."
        )
    elif outlet_pressure is None and inlet_flow is not None:
        warn("No outlet pressure defined with flow inlet. Setting outlet pressure to 0.")


    # Outlet Pressure
    if outlet_pressure:
        if not (isinstance(outlet_pressure, int) or isinstance(outlet_pressure, float)):
            raise TypeError(f"Outlet pressure type '{type(outlet_pressure)}' is not valid. Valid types include float or int.")
        elif not outlet_pressure >= 0:
            raise ValueError(f"Outlet pressure value of '{outlet_pressure}' is not valid. Must be greater than or equal to 0.")

        bcs["outlet"] = {"pressure": outlet_pressure}

    # TODO Add outlet flows potentially?

    # Inlet Pressure
    if inlet_pressure:
        if not (isinstance(inlet_pressure, int) or isinstance(inlet_pressure, float) or isinstance(inlet_pressure, dict)):
            raise TypeError(
                f"Inlet pressure type '{type(inlet_pressure)}' is not valid. Valid types include float or int for single inputs and dict for multiple inputs."
            )
        elif isinstance(inlet_pressure, int) or isinstance(inlet_pressure, float):
            if not inlet_pressure > 0:
                raise ValueError(f"Invalid inlet pressure of '{inlet_pressure}' Pa. Must be greater than 0.")
            bcs["inlet"] = {"pressure": inlet_pressure}
        elif isinstance(inlet_pressure, dict):
            inlet_pressures = {}
            for key, value in inlet_pressure.items():
                if not isinstance(key, int) or not key > 0:
                    raise ValueError(f"Invalid Node Id {key} for multiple inlet pressures. Must be a positive integer.")
                if not (isinstance(value, float) or isinstance(value, int)) or not value > 0:
                    raise ValueError(f"Invalid inlet pressure for node {key} of '{inlet_pressure}' Pa. Must be greater than 0.")
                # Fix indexing here (1-based for IPNODE and IPELEM to 0 based for NetworkX)
                # TODO: Update this if we change input file types.
                inlet_pressures[key - 1] = value
            bcs["inlet"] = {"pressure": inlet_pressures}

    # Inlet flow
    if inlet_flow:
        if not (isinstance(inlet_flow, int) or isinstance(inlet_flow, float) or isinstance(inlet_flow, dict)):
            raise TypeError(
                f"Inlet flow type '{type(inlet_flow)}' is not valid. Valid types include float or int for single inputs and dict for multiple inputs."
            )
        elif isinstance(inlet_flow, int) or isinstance(inlet_flow, float):
            if not inlet_flow > 0:
                raise ValueError(f"Invalid inlet flow of '{inlet_flow}' Pa. Must be greater than 0.")
            bcs["inlet"] = {"flow": inlet_flow}
        elif isinstance(inlet_flow, dict):
            inlet_flows = {}
            for key, value in inlet_flow.items():
                if not isinstance(key, int) or not key > 0:
                    raise ValueError(f"Invalid Element Id {key} for multiple inlet flows. Must be a positive integer.")
                if not (isinstance(value, float) or isinstance(value, int)) or not value > 0:
                    raise ValueError(f"Invalid inlet flow for element {key} of '{inlet_flow}'. Must be greater than 0.")
                # Fix indexing here (1-based for IPNODE and IPELEM to 0 based for NetworkX)
                # TODO: Update this if we change input file types.
                inlet_flows[key - 1] = value
            bcs["inlet"] = {"flow": inlet_flows}

    return bcs


"""
# TODO: Some form that allows:
- inlet type
- inlet BC's

- multiple inlet BC's either flow or pressure (for now say that can only do same type cause otherwise this is hell)
- if multiple inlet BC's:
    - must provide the inlet ID's they are referring to. MATCH THESE ID's with Placentagen.
    NOTE: THIS MEANS WE HAVE TO UPDATE INDEXING SLIGHTLY TO MATCH, AS IPELEM FILES HAVE 1-BASED INDEXING AND NETWORKX HAS 0 BASED INDEXING



- outlet pressure.
- defaults to setting this outlet pressure at all outlet nodes (this makes the most sense to me, can be changed later)



function inputs:

inlet_pressure - value or dict[node_id, pressure]
inlet_flow - value or dict [edge_id, flow]
outlet_pressure - value



function output:

bcs: dict[inlet/outlet, dict[pressure/flow, value or dict[node/element_id, pressure/flow]]]
(if pressure must be node in inner dict, if flow must be element.)

"""
