from fetoflow import *

def main():

    # read in node and element files
    node_filename= 'sample_geometry/FullTree.ipnode'
    element_filename = 'sample_geometry/FullTree.ipelem'

    # process nodes and elements using functions from file_parsing_utils()
    nodes = read_nodes(node_filename)
    elements = read_elements(element_filename)

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
    G = create_geometry(nodes, elements, umbilical_artery_radius, decay_factor, arteries_only, umbilical_vein_radius, decay_factor_vein)
    G = calculate_resistance(G, viscosity_model=viscosity_type)

    # if wanting to incorporate bifurcation pressures losses
    bifurcation_pressure_losses = True

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
    with open(f"output_data/{__file__}_pressure.txt", "w") as f:
        f.write(p)
    with open(f"output_data/{__file__}_flow.txt", "w") as f:
        f.write(q)
if __name__ == '__main__':
    main()
