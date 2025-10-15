import os
import pandas as pd
import numpy as np
import warnings

from .file_parsing_utils import read_nodes, read_elements, define_fields_from_files
from .bc_utils import generate_boundary_conditions
from .geometry_utils import create_geometry, update_geometry_with_pressures_and_flows, create_anastomosis
from .resistance_utils import calculate_resistance
from .matrix_builder import create_matrices,create_small_matrices
from .solve_utils import solve_system,solve_small_system, solve_iterative_system, iterative_solve_small

import os
os.getcwd()


def pressures_and_flows(
    node_filename,
    element_filename,
    boundary_conditions,
    inlet_radius,
    strahler_ratio_arteries,
    *,
    input_directory=".",
    output_directory="./output_data",
    flow_output_filename="flow_values.csv",
    pressure_output_filename="pressure_values.csv",
    arteries_only=False,
    viscosity_model="constant",
    vessel_type="rigid",
    outlet_vein_radius=None,
    strahler_ratio_veins=None,
    anastomosis=None,
    mu=0.33600e-2,  # This is the non-capillary viscosity value used
    capillary_model="analytical2015",
    capillary_parameters=None,  # TODO: Ensure proper dict formatting here.
    radius_filename=None,
    other_field_filenames=None,  # TODO: these are implementable. and radius. Dictionary if used
    verbose=False,
    time_statistics=False,
    return_graph=False,
    # TODO: other input arguments, e.g. mu as a viscosity value? Or capillary settings? All the values of physical stuff and num_generations num_convolutes etc
):
    # Set using laplacian to true and update if needed TODO: Add branching angles here
    use_laplacian = True
    nonlinear = False
    if viscosity_model == "flow_dependent":
        use_laplacian = False
        nonlinear = True

    # Read in data
    if not os.path.exists(input_directory):
        os.makedirs(input_directory)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if not (node_filename.endswith(".ipnode") and element_filename.endswith(".ipelem")):
        raise ValueError("Doesn't work with non ipnode and ipelem files! FIX LATER TODO")  # TODO: fix later

    # TODO: Flesh out these functions and make them more robust with errors.
    nodes = read_nodes(input_directory + os.sep + node_filename)
    elements = read_elements(input_directory + os.sep + element_filename)

    # Boundary Conditions
    inlet_pressure = boundary_conditions.get("inlet_pressure")
    inlet_flow = boundary_conditions.get("inlet_flow")
    outlet_pressure = boundary_conditions.get("outlet_pressure")

    if outlet_pressure is None:
        print("No outlet pressure defined, setting to 0 (gives relative pressure throughout the model).")
        outlet_pressure = 0

    bcs = generate_boundary_conditions(inlet_pressure=inlet_pressure, outlet_pressure=outlet_pressure, inlet_flow=inlet_flow)

    # Create the geometry:
    if not (isinstance(inlet_radius, int) or isinstance(inlet_radius, float)):
        raise ValueError(f"Invalid type '{type(inlet_radius)}' for parameter 'inlet_radius'. Valid types are 'int' and 'float'.")
    if not (isinstance(strahler_ratio_arteries, int) or isinstance(strahler_ratio_arteries, float)):
        raise ValueError(
            f"Invalid type '{type(strahler_ratio_arteries)}' for parameter 'strahler_ratio_arteries'. Valid types are 'int' and 'float'."
        )
    
    fields=None
    # Radius input filenames
    if radius_filename:
        fields = define_fields_from_files({"radius":radius_filename})
    if other_field_filenames:
        if not(isinstance(other_field_filenames, dict)):
            raise ValueError(f"Variable other_field_filenames should be a dict, not type {type(other_field_filenames)}.")
        for other_field, other_field_filename in other_field_filenames.items():
            fields.update(define_fields_from_files({other_field:other_field_filename}))

    if arteries_only:
        G = create_geometry(
            nodes=nodes,
            elements=elements,
            inlet_radius=inlet_radius,
            strahler_ratio_arteries=strahler_ratio_arteries,
            arteries_only=arteries_only,
            fields=fields,
        )
    else:
        if not outlet_vein_radius:
            raise ValueError("No value for parameter 'output_vein_radius' inputted, but are using a venous mesh.")
        if not strahler_ratio_veins:
            raise ValueError("No value for parameter 'strahler_ratio_veins' inputted, but are using a venous mesh.")
        if not (isinstance(outlet_vein_radius, int) or isinstance(outlet_vein_radius, float)):
            raise ValueError(f"Invalid type '{type(outlet_vein_radius)}' for parameter 'outlet_vein_radius'. Valid types are 'int' and 'float'.")
        if not (isinstance(strahler_ratio_veins, int) or isinstance(strahler_ratio_veins, float)):
            raise ValueError(f"Invalid type '{type(strahler_ratio_veins)}' for parameter 'strahler_ratio_veins'. Valid types are 'int' and 'float'.")

        G = create_geometry(
            nodes=nodes,
            elements=elements,
            inlet_radius=inlet_radius,
            strahler_ratio_arteries=strahler_ratio_arteries,
            arteries_only=arteries_only,
            outlet_vein_radius=outlet_vein_radius,
            strahler_ratio_veins=strahler_ratio_veins,
            fields=fields,
        )

    # Create an anastomosis if required: TODO remember IPNODE INDEXING WHEN DOING THE INPUT GRAPH!!
    if anastomosis:
        if not isinstance(anastomosis, dict):
            raise TypeError(f"Variable anastomosis is incorrect type {type(anastomosis)}, should be type(dict).")
        node_from = anastomosis.get("node_from")
        if node_from is None:
            raise ValueError(f"No value for parameter 'node_from' in anastomosis: {anastomosis}")
        node_to = anastomosis.get("node_to")
        if node_to is None:
            raise ValueError(f"No value for parameter 'node_to' in anastomosis: {anastomosis}")
        radius = anastomosis.get("radius")
        # In create_anastomosis, we check the nodes are in the graph and use a dummy radius if one is not inputted. Also converts radius from mm to m.
        G = create_anastomosis(G=G, node_from=node_from, node_to=node_to, radius=radius, mu=mu)

    # Calculate Resistances

    valid_capillary_models = ["analytical2015"]

    if capillary_model not in valid_capillary_models:
        raise ValueError(f"Invalid capillary model '{capillary_model}'. Valid models are {valid_capillary_models}")

    # set up capillary parameters defaults
    num_series = 3
    num_parallel = 6
    num_generations = 3
    num_convolutes_per_branch = 10
    capillary_hematocrit = 0.45
    total_capillary_length = 3.0 / 1000  # m
    capillary_convolute_radius = 7.2e-6  # m
    radius_cap_artery = 1.5e-5  # Intermediate villous tree smallest radius
    radius_cap_scaling_factor = 2  # radius_cap_vein = radius_cap_artery * radius_cap_scaling_factor
    pries_mu = 0.4e-2
    total_segment_length = 1.5 / 1000  # m

    # Pull parameters from 'capillary_parameters' input, and use them instead of defaults if provided.
    if capillary_parameters is None:
        warnings.warn("No capillary parameters passed in. Using default values.", category=UserWarning, stacklevel=2)
    elif not isinstance(capillary_parameters, dict):
        raise TypeError(f"Parameter 'capillary_parameters' type ({type(capillary_parameters)} is invalid. Valid type is 'dict'.)")
    else:
        num_series = capillary_parameters.get("num_series")
        if num_series is None:
            num_series = 3
            capillary_parameters["num_series"] = num_series
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'num_series'. Using default of {num_series}",
                category=UserWarning,
                stacklevel=2,
            )
        if not isinstance(num_series, int):
            raise TypeError(f"Invalid type ({type(num_series)}) for 'capillary_resistance['num_series']'. Valid type is 'int'.")

        num_parallel = capillary_parameters.get("num_parallel")
        if num_parallel is None:
            num_parallel = 6
            capillary_parameters["num_parallel"] = num_parallel
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'num_parallel'. Using default of {num_parallel}",
                category=UserWarning,
                stacklevel=2,
            )
        if not isinstance(num_parallel, int):
            raise TypeError(f"Invalid type ({type(num_parallel)}) for 'capillary_resistance['num_parallel']'. Valid type is 'int'.")

        num_generations = capillary_parameters.get("num_generations")
        if num_generations is None:
            num_generations = 3
            capillary_parameters["num_generations"] = num_generations
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'num_generations'. Using default of {num_generations}",
                category=UserWarning,
                stacklevel=2,
            )
        if not isinstance(num_generations, int):
            raise TypeError(f"Invalid type ({type(num_generations)}) for 'capillary_resistance['num_generations']'. Valid type is 'int'.")

        num_convolutes_per_branch = capillary_parameters.get("num_convolutes_per_branch")
        if num_convolutes_per_branch is None:
            num_convolutes_per_branch = 10
            capillary_parameters["num_convolutes_per_branch"] = num_convolutes_per_branch
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'num_convolutes_per_branch'. Using default of {num_convolutes_per_branch}",
                category=UserWarning,
                stacklevel=2,
            )
        if not isinstance(num_convolutes_per_branch, int):
            raise TypeError(
                f"Invalid type ({type(num_convolutes_per_branch)}) for 'capillary_resistance['num_convolutes_per_branch']'. Valid type is 'int'."
            )

        capillary_hematocrit = capillary_parameters.get("capillary_hematocrit")
        if capillary_hematocrit is None:
            capillary_hematocrit = 0.45
            capillary_parameters["capillary_hematocrit"] = capillary_hematocrit
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'capillary_hematocrit'. Using default of {capillary_hematocrit}",
                category=UserWarning,
                stacklevel=2,
            )
        if not isinstance(capillary_hematocrit, float):
            raise TypeError(f"Invalid type ({type(capillary_hematocrit)}) for 'capillary_resistance['capillary_hematocrit']'. Valid type is 'float'.")

        total_capillary_length = capillary_parameters.get("total_capillary_length")
        if total_capillary_length is None:
            total_capillary_length = 3.0 / 1000
            capillary_parameters["total_capillary_length"] = total_capillary_length
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'total_capillary_length'. Using default of {total_capillary_length * 1000} mm.",
                category=UserWarning,
                stacklevel=2,
            )
        if not (isinstance(total_capillary_length, int) or isinstance(total_capillary_length, float)):
            raise TypeError(
                f"Invalid type ({type(total_capillary_length)}) for 'capillary_resistance['total_capillary_length']'. Valid types are 'int' and 'float'."
            )

        capillary_convolute_radius = capillary_parameters.get("capillary_convolute_radius")
        if capillary_convolute_radius is None:
            capillary_convolute_radius = 7.2e-6  # m
            capillary_parameters["capillary_convolute_radius"] = capillary_convolute_radius
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'capillary_convolute_radius'. Using default of {capillary_convolute_radius * 1000} mm.",
                category=UserWarning,
                stacklevel=2,
            )
        if not (isinstance(capillary_convolute_radius, int) or isinstance(capillary_convolute_radius, float)):
            raise TypeError(
                f"Invalid type ({type(capillary_convolute_radius)}) for 'capillary_resistance['capillary_convolute_radius']'. Valid types are 'int' and 'float'."
            )

        radius_cap_artery = capillary_parameters.get("radius_cap_artery")
        if radius_cap_artery is None:
            radius_cap_artery = 1.5e-5  # Intermediate villous tree smallest radius
            capillary_parameters["radius_cap_artery"] = radius_cap_artery
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'radius_cap_artery'. Using default of {radius_cap_artery * 1000} mm.",
                category=UserWarning,
                stacklevel=2,
            )
        if not (isinstance(radius_cap_artery, int) or isinstance(radius_cap_artery, float)):
            raise TypeError(
                f"Invalid type ({type(radius_cap_artery)}) for 'capillary_resistance['radius_cap_artery']'. Valid types are 'int' and 'float'."
            )

        radius_cap_scaling_factor = capillary_parameters.get("radius_cap_scaling_factor")
        if radius_cap_scaling_factor is None:
            radius_cap_scaling_factor = 2  # radius_cap_vein = radius_cap_artery * radius_cap_scaling_factor
            capillary_parameters["radius_cap_scaling_factor"] = radius_cap_scaling_factor
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'radius_cap_scaling_factor'. Using default of {radius_cap_scaling_factor * 1000} mm.",
                category=UserWarning,
                stacklevel=2,
            )
        if not (isinstance(radius_cap_scaling_factor, int) or isinstance(radius_cap_scaling_factor, float)):
            raise TypeError(
                f"Invalid type ({type(radius_cap_scaling_factor)}) for 'capillary_resistance['radius_cap_scaling_factor']'. Valid types are 'int' and 'float'."
            )

        pries_mu = capillary_parameters.get("pries_mu")
        if pries_mu is None:
            pries_mu = 0.4e-2
            capillary_parameters["pries_mu"] = pries_mu
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'pries_mu'. Using default of {pries_mu} Pa.s.",
                category=UserWarning,
                stacklevel=2,
            )
        if not (isinstance(pries_mu, int) or isinstance(pries_mu, float)):
            raise TypeError(f"Invalid type ({type(pries_mu)}) for 'capillary_resistance['pries_mu']'. Valid types are 'int' and 'float'.")

        total_segment_length = capillary_parameters.get("total_segment_length")
        if total_segment_length is None:
            total_segment_length = 1.5 / 1000
            capillary_parameters["total_segment_length"] = total_segment_length
            warnings.warn(
                f"Parameter 'capillary_paramaters' does not contain a value for 'total_segment_length'. Using default of {total_segment_length * 1000} mm.",
                category=UserWarning,
                stacklevel=2,
            )
        if not (isinstance(total_segment_length, int) or isinstance(total_segment_length, float)):
            raise TypeError(
                f"Invalid type ({type(total_segment_length)}) for 'capillary_resistance['total_segment_length']'. Valid types are 'int' and 'float'."
            )

        # update parameters 'capillary_length' and 'segment_length' as they technically depend on other parameters
        capillary_parameters["capillary_length"] = total_capillary_length / num_parallel
        capillary_parameters["segment_length"] = total_segment_length / num_convolutes_per_branch

    # Calculate resistance
    G = calculate_resistance(G=G, viscosity_model=viscosity_model, mu=mu, capillary_model=capillary_model, capillary_parameters=capillary_parameters)

    # TODO: Potentially record resistance data here if we want it, alongside other data such as radii and length etc?
    ## Good idea as also need resistances to solve small system
    elemData = []
    for u, v in G.edges():
        elemData.append(G.get_edge_data(u,v))
    edge_data_df = pd.DataFrame(elemData)
    edge_data_df = edge_data_df.sort_values(by='edge_id')
    edge_data_df.to_csv(output_directory + os.sep + 'edge_data.csv',index=False)
    # Create A and b matrices:
    if use_laplacian:
            
        ## really dont want to do the svd for a big matrix - it will literally take longer than converging the matrix itself. im trying to find a good way
        ## to approximate it for my gmres switching.

        A, b, bc_export, _ = create_small_matrices(G=G, bcs=bcs)
    else:
        A, b = create_matrices(G=G, n=G.number_of_nodes(), m=G.number_of_edges(), bcs=bcs)
        
    # Solve Matrix
    if not use_laplacian:
        if nonlinear:
            pressures, flows = solve_iterative_system(G, A, b, G.number_of_nodes(), G.number_of_edges(), bcs, viscosity_model, mu, capillary_model, capillary_parameters)
        else: # Linear
            pressures, flows = solve_system(A=A, b=b, num_nodes=G.number_of_nodes(), num_edges=G.number_of_edges())
    else:
        if nonlinear:
            pressures, flows = iterative_solve_small(A=A, b=b, G=G, boundary_conditions=bc_export)
        else:
            pressures, flows, _ = solve_small_system(A=A, b=b, G=G, boundary_conditions=bc_export)


    # Export to output directories. TODO: All the other stuff we want to export. Will be more done in the solver function.
    pressures = pd.DataFrame([{"node": i, "pressure": pressures[i]} for i in G.nodes()])
    pressures.to_csv(output_directory + os.sep + pressure_output_filename, index=False)
    flows = pd.DataFrame([{"element": data["edge_id"], "flow": flows[data["edge_id"]]} for u, v, data in G.edges(data=True)])
    flows = flows.sort_values(by="element")
    flows.to_csv(output_directory + os.sep + flow_output_filename, index=False)

    # Add pressure and flow fields to the newtorkX graph object, which we can return if wanted for a more in-depth investigation
    if return_graph:
        G = update_geometry_with_pressures_and_flows(G, pressures, flows)
        return G
    return
