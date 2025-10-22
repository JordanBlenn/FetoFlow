from FetoFlow import *
import pandas as pd

def main():

    # read in node and element files
    node_filename= 'sample_geometry/P49_full_tree.ipnode'
    element_filename = 'sample_geometry/P49_full_tree.ipelem'
    radius_filename = 'sample_geometry/P49_full_tree_radius.ipfiel'

    # process nodes and elements using functions from file_parsing_utils()
    nodes = read_nodes(node_filename)
    elements = read_elements(element_filename)
    fields = define_fields_from_files(files={"radius":radius_filename})

    # define boundary conditions
    inlet_pressure, outlet_pressure = 6650, 2660

    # generate boundaray conditions dictionary 
    bcs = generate_boundary_conditions(inlet_pressure = inlet_pressure, outlet_pressure = outlet_pressure, inlet_flow=None)

    # define other required geometric features (radii and decay factors)
    umbilical_artery_radius, decay_factor = 1.8 / 1000, 1.38 
    umbilical_vein_radius, decay_factor_vein = 4.0 / 1000, 1.46
    arteries_only = False # this should rarely be true

    viscosity_type = 'constant' # can also be 'pries_network' or 'pries_vessel' if wanting to incorporate radius-dependence

    # Generate the di-graph & calculate the resistances based on the viscosity
    G = create_geometry(nodes, elements, umbilical_artery_radius, decay_factor, arteries_only, umbilical_vein_radius, decay_factor_vein,fields=fields)
    radii_out = getRadii(G)
    radii_export = pd.DataFrame([{"Node In" : node_in, "Node Out" : node_out, "Radius" : r} for (node_in,node_out),r in radii_out.items()])
    radii_export.to_csv("output_data/radius_from_file.csv")
    G = calculate_resistance(G, viscosity_model=viscosity_type)

    # if wanting to incorporate bifurcation pressures losses
    bifurcation_pressure_losses = False

    # flow-dependent visc?????

    if bifurcation_pressure_losses:
        # calculate branch angles for use in bifurcation pressure losses
        calculate_branching_angles(G)
        # create small matrices, knowing that we need to solve iteratively
        A,b,bc_export,iterative_info = create_small_matrices(G,bcs,branching_angles=True)

        # define convergence tolerance and solve system iteratively
        convergence_tol = 0.01
        p,q = iterative_solve_small(A,b,G,bc_export,convergence_tol,iterative_info,use_gmres=False,adaptive_stepping=True)
    else:
        # solve directly
         A,b,bc_export = create_small_matrices(G,bcs,branching_angles=False)
         p,q = solve_small_system(A,b,G,bc_export)

    # store pressure and flows
    pressures = pd.DataFrame([{"Node" : node, "Pressure" : pressure} for node,pressure in p.items()])
    pressures.to_csv("output_data/example_simulation_pressures.csv")
    flows = pd.DataFrame([{"Element" : element, "Flow" : flow} for element,flow in q.items()])
    flows.to_csv("output_data/example_simulation_flows.csv")

if __name__ == '__main__':
    main()
