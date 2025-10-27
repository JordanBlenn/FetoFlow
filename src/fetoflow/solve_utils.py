import scipy.sparse as sps
import numpy as np
import time
from scipy.special import expit

from .matrix_builder import create_matrices
from .resistance_utils import calculate_resistance, calculate_viscosity_factor_from_radius
from .geometry_utils import update_geometry_with_pressures_and_flows

def solve_small_system(A, b, G, boundary_conditions, ill_conditioned=False):
    """Solve a reduced system of equations for internal pressures and flows.

    This function solves the reduced system obtained after eliminating known boundary 
    values. For flow boundary conditions, it reconstructs inlet pressures from the
    solution.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        Reduced system matrix after boundary elimination
    b : numpy.ndarray
        Right-hand side vector after boundary elimination
    G : networkx.DiGraph
        Graph containing network topology and edge properties
    boundary_conditions : tuple
        (bc_type, boundary_indices, boundary_vals, inlet_idx) describing boundary
        conditions and indices
    ill_conditioned : bool, optional
        If True, return internal pressures for iterative refinement
    p0 : numpy.ndarray, optional
        Initial pressure guess for iterative methods
    current_p : numpy.ndarray, optional 
        Current pressure solution for warm starts
    max_iterations : int, optional
        Maximum iterations for iterative solver
    restart : int, optional
        GMRES restart parameter

    Returns
    -------
    dict
        Mapping of node indices to pressures
    dict
        Mapping of edge indices to flows
    numpy.ndarray, optional
        Internal pressure solution if ill_conditioned=True
    """
    # Unpack boundary conditions tuple
    bc_type, boundary_indices, boundary_vals, inlet_idx = boundary_conditions

    # Direct solve of reduced system Ap = b for internal pressures
    p = sps.linalg.spsolve(A, b)  
    
    # Save internal pressures if needed for iterative refinement
    if ill_conditioned:
        internal_p = p.copy()
    
    # Initialize flow array for all edges
    q = np.zeros(shape=G.number_of_edges())

    # For flow BCs, reconstruct inlet pressures using p = QR + p_downstream
    if bc_type == "Flow":
        boundary_vals_current_iteration = boundary_vals.copy()  # Make copy to avoid modifying original
        for current_inlet in inlet_idx:
            # Get downstream node and flow value
            adj_to_inlet = list(G.out_edges(current_inlet))[0][1]
            value_idx = np.where(boundary_indices == current_inlet)[0][0]
            # Account for removed inlet nodes in index mapping
            index_adjustment = np.sum([adj_to_inlet > i for i in inlet_idx])
            # Reconstruct inlet pressure: p_inlet = Q*R + p_downstream
            p0 = boundary_vals[value_idx]*G[current_inlet][adj_to_inlet]["resistance"] + p[adj_to_inlet - index_adjustment]
            boundary_vals_current_iteration[value_idx] = p0 
    else:
        # For pressure BCs, use given boundary values directly
        boundary_vals_current_iteration = boundary_vals

    # Insert boundary values back into pressure vector in correct order
    indices = np.argsort(boundary_indices)
    for idx,val in zip(boundary_indices[indices],boundary_vals_current_iteration[indices]):
        p = np.insert(p,idx,val)
    # Calculate flows using Q = (p1-p2)/R for each edge
    for u,v in G.edges():
        pu,pv = p[u],p[v]  # Get pressures at each end
        q[G[u][v]['edge_id']] = (pu-pv)/G[u][v]['resistance']  # Flow from pressure difference

    # Pack results into dictionaries with node/edge indices as keys
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    pressures = {node_id: p[node_id] for node_id in range(num_nodes)}
    flows = {edge_id: q[edge_id] for edge_id in range(num_edges)}

    # Return appropriate results based on ill_conditioned flag
    if ill_conditioned:
        return pressures, flows, internal_p  # Include internal pressures for refinement
    return pressures, flows  # Standard return of results

def __solve_with_gmres(A, b, G, boundary_conditions, current_p=None):
    """Solve system using GMRES with ILU preconditioning.
    
    Internal helper that implements GMRES solution with incomplete LU 
    preconditioning for ill-conditioned systems. Uses scipy.sparse.linalg.gmres
    with an ILU preconditioner.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        System matrix 
    b : numpy.ndarray
        Right-hand side vector
    G : networkx.DiGraph
        Graph containing network topology
    boundary_conditions : tuple
        (bc_type, boundary_indices, boundary_vals, inlet_idx)
    current_p : numpy.ndarray, optional
        Initial guess for GMRES

    Returns
    -------
    dict
        Node pressures
    dict 
        Edge flows
    numpy.ndarray
        Internal pressure solution
    int
        GMRES convergence flag (0 = success)
    """
    bc_type, boundary_indices, boundary_vals, inlet_idx = boundary_conditions
    # incomplete lu conditioner as per Al-Kurdi/Kincaid 
    n = A.shape[0]
    A_inv = sps.linalg.spilu(A=A,drop_tol=np.min(A.diagonal()),fill_factor=15) # approximates A inverse using an incomplete LU factorisation
    M = sps.linalg.LinearOperator(shape=(n,n),matvec=A_inv.solve)
    v = np.random.rand(n)
    print(f"Approximate Preconditioner Performance ||MAx - x|| for random x: {np.linalg.norm(M @ (A @ v) - v)}")
    residuals = []
    def callback(rk):
        residuals.append(rk)
    b = b.reshape(b.shape[0]) # look at this!
    p,flag = sps.linalg.gmres(A=A,b=b,x0=current_p,M=M,maxiter=300,callback=callback,restart=100)
    if flag != 0:
        print("gmres not converged")
        solver_residual = residuals[-1]
        print(f"current (avg) residual: {solver_residual/len(b)}")
    else:
        print("gmres converged")
    internal_p = p.copy()
    q = np.zeros(shape=G.number_of_edges())
    if bc_type == "Flow":
        for current_inlet in inlet_idx:
            adj_to_inlet = list(G.out_edges(current_inlet))[0][1]
            value_idx = np.where(boundary_indices == current_inlet)[0][0]
            p0 = boundary_vals[value_idx]*G[current_inlet][adj_to_inlet]["resistance"] + p[current_inlet]
            boundary_vals[value_idx] = p0

    indices = np.argsort(boundary_indices)
    for idx,val in zip(boundary_indices[indices],boundary_vals[indices]):
        p = np.insert(p,idx,val)
    for u,v in G.edges():
        pu,pv = p[u],p[v]
        q[G[u][v]['edge_id']] = (pu-pv)/G[u][v]['resistance']
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    pressures = {node_id: p[node_id] for node_id in range(num_nodes)}
    flows = {edge_id: q[edge_id] for edge_id in range(num_edges)} #TODO: triple check these but they should be fine
    return pressures, flows,internal_p,flag

    
def update_small_matrix(G, diag_to_update, flows, iter_options):
    """Update matrix for bifurcation pressure losses.

    Updates the diagonal entries of the system matrix to account for pressure 
    losses at vessel bifurcations. Handles both diverging (1-to-many) and 
    converging (many-to-1) bifurcations.

    Parameters
    ----------
    G : networkx.DiGraph
        Graph containing network topology and vessel properties 
    diag_to_update : numpy.ndarray
        Diagonal entries of the matrix to update with loss terms
    flows : dict
        Edge flows from previous iteration
    iter_options : dict
        Dictionary containing options for iteration:
        - branch_nodes: Mapping of branch node indices to (entering_nodes, leaving_nodes)
        - inlet_bc: Boundary condition type ("Flow" or "Pressure")
        - max_inlet_bc_flow: Maximum inlet flow value for flow BCs

    Returns
    -------
    numpy.ndarray
        Updated diagonal entries including bifurcation losses

    Notes
    -----
    Implements the Mynard bifurcation pressure loss model with Reynolds number
    scaling for low Re flows. The loss coefficient depends on:
    - Flow ratio between parent and daughter vessels
    - Area ratio between vessels
    - Bifurcation angles
    - Reynolds number based viscous effects
    """
    # Start timing for performance tracking
    t = time.time()
    # depending on diverging/converging use associated edges to compute pressure drop coefficients for the current bifurcation/edge combo.
    # could look at only computing some of these once for efficiency, but i think since the flows need to be recalculated and updated anyway not worth
    t = time.time()  # Start timing for performance tracking
    
    if iter_options["branch_nodes"]:
        # Thresholds for Reynolds number scaling
        viscous_re_threshold = 100   # Below this, scale losses by Re/threshold
        absolute_re_threshold = 300   # Max allowed Re (umbilical vessels limit)
        
        # Blood properties
        rho = 1060      # Density (kg/m^3)
        mu = 3.36e-3    # Dynamic viscosity (Pa.s)
        branching_dict = iter_options["branch_nodes"]
        if iter_options["inlet_bc"] == "Pressure":
            flow_inlet = np.inf
        else:
            flow_inlet = iter_options["max_inlet_bc_flow"]
        for row in branching_dict.keys():
            # get required info
            entering_nodes,leaving_nodes,original_row = branching_dict[row]
            if len(leaving_nodes) > 1:
                datum_flow,datum_area = flows[G[entering_nodes[0]][original_row]["edge_id"]],np.pi*G[entering_nodes[0]][original_row]["radius"]**2
                max_flow = min(absolute_re_threshold*datum_area*mu/(rho*2*np.pi*G[entering_nodes[0]][original_row]["radius"]),flow_inlet)
                datum_flow = np.clip(datum_flow,1e-14,max_flow)
                for j in leaving_nodes:
                    flow_j,area_j = flows[G[original_row][j]["edge_id"]],np.pi*G[original_row][j]["radius"]**2
                    max_flow = min(absolute_re_threshold*area_j*mu/(rho*2*np.pi*G[original_row][j]["radius"]),flow_inlet)
                    flow_j = np.clip(flow_j,1e-14,max_flow)

                    coef = 1 - 1/((flow_j/datum_flow)*(datum_area/area_j))*np.cos((3/4*np.pi - 3/4*(np.pi - G[original_row][j]["bifurcation_angle"])))
                    p_loss = coef*rho*(flow_j/area_j)**2 + 1/2*rho*((datum_flow/datum_area)**2 - (flow_j/area_j)**2) # full pressure loss as per Mynard

                    Re = rho * (flow_j/area_j)* 2 * G[original_row][j]["radius"]/mu
                    if Re < viscous_re_threshold:
                        p_loss*= Re/viscous_re_threshold
                    effective_res = p_loss/flow_j
                    if effective_res > 1/diag_to_update[G[original_row][j]["edge_id"]]:
                        # print(f"Predicted unphysical pressure loss in for leaving edge {G[original_row][j]["edge_id"]}, loss has been constrained to increase resistance by 20%")
                        effective_res = 0.2/diag_to_update[G[original_row][j]["edge_id"]] # this is only really a contingency for flow boundary conditions where 
                    diag_to_update[G[original_row][j]["edge_id"]] = 1/(1/diag_to_update[G[original_row][j]["edge_id"]] + effective_res)

            elif len(entering_nodes) > 1:
                datum_flow,datum_area = flows[G[original_row][leaving_nodes[0]]["edge_id"]],np.pi*G[original_row][leaving_nodes[0]]["radius"]**2
                max_flow = min(absolute_re_threshold*datum_area*mu/(rho*2*np.pi*G[entering_nodes[0]][original_row]["radius"]),flow_inlet)
                datum_flow = np.clip(datum_flow,1e-14,max_flow)
                # Initialize total pressure loss for converging bifurcation
                p_loss = 0
                for j in entering_nodes:
                    # Get flow and cross-sectional area for entering vessel
                    flow_j = flows[G[j][original_row]["edge_id"]]
                    area_j = np.pi*G[j][original_row]["radius"]**2
                    
                    # Limit flow to prevent unphysical Reynolds numbers
                    max_flow = min(absolute_re_threshold*area_j*mu/(rho*2*np.pi*G[j][original_row]["radius"]), flow_inlet)
                    flow_j = np.clip(flow_j, 1e-14, max_flow)
                    
                    # Compute loss coefficient based on Mynard model
                    coef = 1 - 1/((flow_j/datum_flow)*(datum_area/area_j))*np.cos((3/4*np.pi - 3/4*(np.pi - G[j][original_row]["bifurcation_angle"])))
                    
                    # Add contribution to total pressure loss
                    p_loss += coef*rho*(flow_j/area_j)**2 + 1/2*rho*((datum_flow/datum_area)**2 - (flow_j/area_j)**2)
                    
                    # Scale losses by Reynolds number if viscous forces dominate
                    Re = rho * (flow_j/area_j)* 2* G[j][original_row]["radius"]/mu
                    if Re < viscous_re_threshold:
                        p_loss *= Re/viscous_re_threshold  # Linear scaling with Re
                # Convert total pressure loss to effective resistance
                effective_res = p_loss/datum_flow
                
                # Prevent unphysical resistance increases
                if effective_res > 1/diag_to_update[G[original_row][leaving_nodes[0]]["edge_id"]]:
                    # print(f"Predicted unphysical pressure loss for leaving edge {G[original_row][leaving_nodes[0]]["edge_id"]}, loss has been constrained to increase resistance by 20%")
                    effective_res = 0.2/diag_to_update[G[original_row][leaving_nodes[0]]["edge_id"]]
                
                # Update resistance 
                diag_to_update[G[original_row][leaving_nodes[0]]["edge_id"]] = 1/(1/diag_to_update[G[original_row][leaving_nodes[0]]["edge_id"]] + effective_res)
    
    # Report timing statistics
    print(f"matrix update: {time.time() - t}")
    return diag_to_update




def solve_system(A, b, num_nodes, num_edges):
    """Solve full system of equations for pressures and flows.
    
    Solves the full system obtained from create_matrices() that includes 
    both pressure and flow equations. The solution vector contains both
    nodal pressures and edge flows.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        Full system matrix with pressure and flow equations
    b : numpy.ndarray
        Right-hand side vector
    num_nodes : int
        Number of nodes in network
    num_edges : int
        Number of edges in network

    Returns
    -------
    dict
        Mapping of node indices to pressures
    dict
        Mapping of edge indices to flows

    Notes
    -----
    The solution vector x contains:
    - First num_nodes entries: nodal pressures
    - Last num_edges entries: edge flows
    """
    # Direct solve of full system using sparse solver
    x = sps.linalg.spsolve(A, b)
    
    # Extract pressures from first section of solution vector
    pressures = {node_id: x[node_id] for node_id in range(num_nodes)}
    
    # Extract flows from second section, offset by number of nodes
    flows = {edge_id: x[num_nodes + edge_id] for edge_id in range(num_edges)}
    
    return pressures, flows

def iterative_solve_small(A, b, G, bc_export, tol, info, alpha=1, maxiter=20, use_gmres=False, adaptive_stepping=True):
    """Iteratively solve system with bifurcation pressure losses.

    Solves the system iteratively to handle nonlinear bifurcation pressure
    losses. Uses either direct solution or GMRES with optional adaptive stepping
    for better convergence.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix  
        Initial system matrix
    b : numpy.ndarray
        Right-hand side vector
    G : networkx.DiGraph
        Graph containing network topology
    bc_export : tuple
        (bc_type, boundary_indices, boundary_vals, inlet_idx)
    tol : float
        Convergence tolerance on pressure MSE
    info : dict
        Options dictionary containing:
        - branching_update_matrix: Initial resistance matrix
        - branching_calc_matrices: Pre-computed matrices for efficiency
        - inlet_bc: Boundary condition type 
        - max_inlet_bc_flow: Max inlet flow for flow BCs
    alpha : float, optional
        Initial relaxation factor
    maxiter : int, optional  
        Maximum iterations
    use_gmres : bool, optional
        Whether to use GMRES solver
    adaptive_stepping : bool, optional
        Whether to adapt alpha based on convergence
        
    Returns
    -------
    dict
        Node pressures
    dict
        Edge flows

    Notes
    -----
    Uses a relaxed fixed-point iteration:
    1. Solve system with current resistances
    2. Update bifurcation loss terms
    3. Relaxed update of solution
    4. Check convergence
    """
    # Start timing
    s = time.time()
    
    # Get initial resistance matrix
    W = info["branching_update_matrix"]
    
    # Initial solve using appropriate method
    if use_gmres:
        p, q, internal_p = solve_small_system(A, b, G, bc_export, ill_conditioned=True)
    else:
        p, q = solve_small_system(A, b, G, bc_export, ill_conditioned=False)
    
    # Convert dictionary results to arrays for computation
    p0 = np.array([p[node] for node in p.keys()])
    q0 = np.array([q[elem] for elem in q.keys()])
    
    # Print initial solution statistics
    print(f"Info |  Max Pressure: {np.max(p0)} | Min Pressure: {np.min(p0)} | Max flow: {np.max(q0)} | Min flow: {np.min(q0)}")

    # Initialize convergence tracking
    p_mse = np.inf  # Mean squared error in pressures
    old_mse = np.inf  # Previous iteration MSE
    iteration = 1
    
    print("Warm up period...")
    # Set relaxation factor bounds
    min_alpha = min(5*tol, 1)  # Lower bound proportional to tolerance
    max_alpha = 1              # Upper bound at full update
    flag = int(use_gmres)     # Track GMRES convergence
    init_alpha = alpha        # Store initial alpha for adaptive stepping
    
    # Extract pre-computed matrices based on problem type
    if len(info["branching_calc_matrices"]) == 1:
        # Simple case - just need reduced incidence matrix
        Br = info["branching_calc_matrices"][0]
    else:
        # Full flow BC case - need additional matrices
        [Br, u, c, qi, vp_out, edge_idx] = info["branching_calc_matrices"]
    
    # Main iteration loop
    while p_mse > tol or flag:
        # Check iteration limit
        if iteration > maxiter:
            print("Maximum Iteration Count Reached. Consider constraining the learning rate parameter alpha (if observing cycling) or GMRES if MSE is very large")
            break
            
        # Update resistances based on current flows
        diag_update = update_small_matrix(G, W.diagonal(), q, info)
        W_new = sps.diags(diag_update, offsets=0).tocsc()
        
        # Reconstruct system matrix based on problem type
        if len(info["branching_calc_matrices"]) == 1:
            # pressure BC
            A = Br @ W_new @ Br.T
        else:
            # Flow BC
            a_inv = sps.diags(diag_update[edge_idx]).tocsc()
            A = Br @ W_new @ Br.T - u @ a_inv @ c
            b = -u @ a_inv @ qi - vp_out
        if use_gmres:
            p,q,internal_p,flag = __solve_with_gmres(A,b,G,bc_export,current_p=internal_p)
        else:
            p,q = solve_small_system(A,b,G,bc_export,ill_conditioned=False)
        p_init = np.array([p[node] for node in p.keys()])
        p_mse = 1/len(p_init)*np.linalg.norm(p_init - p0)**2

        if adaptive_stepping:
            if old_mse != np.inf and iteration > maxiter//10: # let the system warmup for a sec (10% of total iter)
                if p_mse < old_mse: # converging
                    alpha = np.clip(p_mse/old_mse,min_alpha,max_alpha)
                else: # diverging - be more carful!
                    alpha = np.clip(old_mse/p_mse,min_alpha,init_alpha)

        p1 = (1-alpha)*p0 + alpha*p_init # accepted p solution based on alpha
        q1 = (1-alpha)*q0 + alpha*np.array([q[elem] for elem in q.keys()])

        p_diff = p1 - p0
        q_diff = q1 - q0
        p_mse = 1/len(p)*np.linalg.norm(p_diff)**2
        q_mse = 1/len(q)*np.linalg.norm(q_diff)**2
        p_infnorm, q_infnorm = np.max(np.abs(p_diff)),np.max(np.abs(q_diff))
        p0 = p1
        q0 = q1
        print(f"Info |  Max Pressure: {np.max(p0)} | Min Pressure: {np.min(p0)} | Max flow: {np.max(q0)} | Min flow: {np.min(q0)}")
        print(f"Iteration: {iteration} | Pressure MSE:  {round(p_mse,4)} | Flow MSE: {q_mse}, | alpha : {alpha} | Max diff pressure: {p_infnorm} | Max diff flow: {q_infnorm}")
        if iteration == maxiter//10: 
            print("Warm up period COMPLETE")
            if adaptive_stepping:
                print("Adaptive stepping STARTING")
        old_mse = p_mse
        iteration += 1
    print(f"solve time: {time.time() - s}")
    return p,q


def solve_iterative_system(G, A, b, num_nodes, num_edges, bcs, viscosity_model, mu, capillary_model, capillary_parameters, tol=0.01, max_solve_time=120):
    """Solve full system iteratively with nonlinear effects.
    
    Iteratively solves the full system accounting for vessel elasticity, 
    non-linear blood rheology, and geometry updates. Uses direct solution
    with resistance updates until convergence.

    Parameters
    ----------
    G : networkx.DiGraph
        Graph containing network topology
    A : scipy.sparse.csr_matrix
        Initial system matrix
    b : numpy.ndarray 
        Right-hand side vector
    num_nodes : int
        Number of nodes in network
    num_edges : int
        Number of edges in network
    bcs : dict
        Boundary conditions dictionary
    viscosity_model : str
        Type of viscosity model to use
    mu : float
        Base blood viscosity
    capillary_model : str
        Type of capillary model
    capillary_parameters : dict
        Parameters for capillary model
    tol : float, optional
        Convergence tolerance, default 0.01
    max_solve_time : float, optional
        Maximum solve time in seconds, default 120s

    Returns
    -------
    dict
        Node pressures
    dict
        Edge flows

    Notes
    -----
    Algorithm:
    1. Solve system directly
    2. Update geometry and properties
    3. Update resistances
    4. Check convergence
    5. Repeat until converged or timeout

    Current limitations:
    - Does not handle branching angles
    - Basic viscosity models only
    - No coupled wall mechanics
    """
    # TODO: Make solve iterative small system work as well
    # Solve matrix directly, update resistances, check convergence.
    # Used for elasticity, non-linear (flow-dependent) blood rheology, and branching-angle effect.
    import copy
    x = sps.linalg.spsolve(A, b)
    A_old = copy.deepcopy(A)
    b_old = copy.deepcopy(b)

    pressures = {node_id: x[node_id] for node_id in range(num_nodes)}
    flows = {edge_id: x[num_nodes + edge_id] for edge_id in range(num_edges)}
    G = update_geometry_with_pressures_and_flows(G, pressures, flows)
    G = update_graph(G,  viscosity_model, mu, capillary_model, capillary_parameters)
    
    A, b = create_matrices(G, num_nodes, num_edges, bcs)

    # A and A_old are csr_matrix
    # diff_data = np.abs(A.data - A_old.data)
    # max_diff = diff_data.max()
    # print("Max difference in nonzero entries:", max_diff)

    # if max_diff > 1e-12:
    #     print("A has changed")
    # else:
    #     print("A is effectively unchanged")
    # raise ValueError

    # A = update_A_matrix(A, pressures, flows)
    # Convergence - also timeout time limit.
    solve_start_time = time.time()
    current_time = time.time()
    iteration_counter = 0
    while current_time - solve_start_time < max_solve_time:
        iteration_counter += 1
        x_new = sps.linalg.spsolve(A, b)
        if np.linalg.norm(x_new - x) < tol:
            # Converges
            print(f"Non-linear solution converged! {iteration_counter} iterations.")
            break # This uses the old iteration values currently, don't think we need to use the new ones if it converges
        # Pull out pressures and flows as required for updating A Matrix
        pressures = {node_id: x[node_id] for node_id in range(num_nodes)}
        flows = {edge_id: x[num_nodes + edge_id] for edge_id in range(num_edges)}
        # Update G, A and x
        G = update_geometry_with_pressures_and_flows(G, pressures, flows)
        G = update_graph(G, viscosity_model, mu, capillary_model, capillary_parameters)
        A, b = create_matrices(G, num_nodes, num_edges, bcs)
        x = x_new
        current_time = time.time()
        iteration_counter += 1
    # If timed out, print msg.
    if current_time - solve_start_time >= max_solve_time:
        print(f"Solution timed out before convergence (time limit of {max_solve_time} seconds)! Returning values from last iteration...")

    return pressures, flows

def update_A_matrix(A, pressures, flows, G, n, m, flow_dependent_viscosity=False, branching_angles=False, elastic_vessels=False):
    """Update system matrix with nonlinear effects.
    
    DEPRECATED: This function is incomplete and not recommended for use.
    Use the small matrix formulation instead for branching angles and 
    other nonlinear effects.

    Parameters
    ----------
    A : scipy.sparse.csr_matrix
        System matrix to update
    pressures : dict
        Node pressures 
    flows : dict
        Edge flows
    G : networkx.DiGraph
        Network graph
    n : int
        Number of nodes
    m : int
        Number of edges
    flow_dependent_viscosity : bool, optional
        Whether to include flow-dependent viscosity
    branching_angles : bool, optional
        Whether to include branching angle effects
    elastic_vessels : bool, optional
        Whether to include vessel elasticity effects

    Returns
    -------
    scipy.sparse.csr_matrix
        Updated system matrix

    Warnings
    --------
    This function is deprecated and incomplete. Use small matrix 
    formulation in create_small_matrices() instead.
    """
    # TODO: Write this function
    if branching_angles:
        raise Warning("Dont use this for branching angles! Use the small matrix system!")
        rho = 1060 # kg/m^3
        for row in range(m):
            flow_cons_row = A[row,m:m+n]
            if flow_cons_row.nnz == 3: # bifurcation

                leaving_edges = list(G.out_edges(row))
                entering_edges = list(G.in_edges(row))
                leaving_nodes = [leaving_edges[i][1] for i in range(len(leaving_edges))]
                entering_nodes = [entering_edges[i][0] for i in range(len(entering_edges))]

                if len(leaving_nodes) > 1:
                    datum_flow,datum_area = flows[G[entering_nodes[0]][row]["edge_id"]],np.pi*G[entering_nodes[0]][row]["radius"]**2
                    for j in leaving_nodes:
                        flow_j,area_j = flows[G[row][j]["edge_id"]],np.pi*G[row][j]["radius"]**2
                        coef = 1 - 1/((flow_j/datum_flow)*(datum_area/area_j))*np.cos((3/4*np.pi - 3/4*(np.pi - G[row][j]["bifurcation_angle"])))
                        A[m+G[row][j]["edge_id"],row] = 1 - max(coef*1/2*rho*(datum_flow/datum_area)**2/pressures[row],0)
                elif len(entering_nodes) > 1:

                    datum_flow,datum_area = flows[G[row][leaving_nodes[0]]["edge_id"]],np.pi*G[row][leaving_nodes[0]]["radius"]**2
                    coef = 1
                    for j in entering_nodes:
                        flow_j,area_j = flows[G[j][row]["edge_id"]],np.pi*G[j][row]["radius"]**2
                        coef -= 1/((flow_j/datum_flow)*(datum_area/area_j))*np.cos((3/4*np.pi - 3/4*(np.pi - G[j][row]["bifurcation_angle"])))
                    A[m+G[row][leaving_nodes[0]]["edge_id"],row] = 1 - max(coef*1/2*rho*(datum_flow/datum_area)**2/pressures[row],0) # lower pressures -> 0

    return A

def update_graph(G, viscosity_model, mu, capillary_model, capillary_parameters, flow_dependent_viscosity=False):
    """Update graph with hematocrit and viscosity changes.
    
    Updates vessel properties based on the Pries et al. (1990) plasma 
    skimming model and recalculates resistances. Handles hematocrit 
    distribution at bifurcations and corresponding viscosity changes.

    Parameters
    ----------
    G : networkx.DiGraph
        Graph containing network topology and vessel properties
    viscosity_model : str
        Type of viscosity model to use
    mu : float
        Base blood viscosity
    capillary_model : str
        Type of capillary model to use
    capillary_parameters : dict
        Parameters for capillary model
    flow_dependent_viscosity : bool, optional
        Whether to use flow-dependent viscosity

    Returns
    -------
    networkx.DiGraph
        Updated graph with new hematocrit and resistance values

    Notes
    -----
    Algorithm:
    1. Calculate fractional RBC flow (FQE) for each daughter vessel
    2. Update daughter hematocrit based on parent values
    3. Recalculate viscosity factors
    4. Update vessel resistances

    Based on:
    Pries et al. (1990) - Blood viscosity in tube flow: dependence on 
    diameter and hematocrit. Am J Physiol.

    Limitations:
    - Multiple parent vessels not supported
    - Uses mean diameter for >2 daughter vessels
    - Assumes constant plasma viscosity
    """
    # Steps:
    # Define X0, A, B.
    # Work out the FQE. fraction of red blood cells in daughter cell.
    # Then new hematocrit = FQE/FQB * previous hematocrit.
    # Then update resistances based on hematocrit.
    # Then return graph.
    # This stuff based on pries paper 1990

    # THe pries paper is shit, and I'm going to be doing things here that may or may not work but we'll see what happens.
    eps = 1e-9 # need a tiny number for clipping.
    # THIS FUCNKING USES MICROMETERS AS UNITS I HATE IT SO MUCH!

    # Loop through vessels
    for u, v, _ in G.edges(data=True):
        # Get the parent node to see where the flow goes
        parent_node = list(G.predecessors(u))
        if not parent_node:
            continue # Don't update anything here? Cause inlet so flow hasn't broken yet.
        elif len(parent_node) > 1:
            continue #This breaks down with multiple parents, so I'm going to pretend it doesn't exist for now.
        FQ_B =  G[u][v]['flow'] / G[parent_node[0]][u]['flow']

        # Get children nodes
        child_nodes = list(G.successors(u))
        assert len(child_nodes) > 0 # Should always at least have V as a child.
        # Remove v from child nodes
        child_nodes.remove(v)
        if len(child_nodes) == 0: # Bifurcation doesn't properly split into 2, so leave as it is. Hematocrit should be the same as parent node.
            G[u][v]['hematocrit'] =  H_d
            continue

        # Parent hematocrit
        H_d = G[parent_node[0]][u]['hematocrit']
        H_d = np.clip(H_d, eps, 1 - eps) # Prevent Hematocrit from breaching physical limits where approximation breaks down.

        # Diameters:
        D_f = G[parent_node[0]][u]['radius'] * 2 * 10**6 # m to um
        D_alpha = G[u][v]['radius'] * 2 * 10**6 # m to um
        # For D_beta, pries formulation assumes only 2 daughter vessels. 
        # To take this into account, if there are more than 2, we take the average diameter.
        
        D_beta = sum(G[u][child_node]['radius'] * 2 for child_node in child_nodes) / len(child_nodes) * 10 ** 6 # Average of other children, m-um.
        
        X0 = 0.4 / D_f
        assert X0 < 0.5, f"Maths breaks down otherwise, D_f={D_f}"
        A = -6.96*np.log(D_alpha/D_beta)/D_f
        B = 1 + 6.98*(1-H_d)/D_f # Changing D_f into um for these equations cause I think that's what they based them on.

        FQ_E = expit(A + B * (np.log((FQ_B-X0)) - np.log(1 - FQ_B - X0 ))) # Rewritten here for stability.

        G[u][v]['hematocrit'] = np.clip(H_d * FQ_E / FQ_B, eps, 1- eps) # Clip to physical values again
        # if(G[u][v]['hematocrit']) == 0: raise ValueError("BREAKING")
        # print(FQ_E, FQ_B)
        G[u][v]['viscosity_factor'] = calculate_viscosity_factor_from_radius(G[u][v]['radius'], G[u][v]['hematocrit'])

    # Update resistances. TODO Update viscosity model, mu, capillary model and capillary parameters from defaults!
    import copy
    G_copy = copy.deepcopy(G)
    G = calculate_resistance(G=G, viscosity_model=viscosity_model, mu=mu, capillary_model=capillary_model, capillary_parameters=capillary_parameters)
    print("Nodes equal:", G.nodes(data=True) == G_copy.nodes(data=True))
    print("Edges equal:", G.edges(data=True) == G_copy.edges(data=True))
    # differences = []

    # for u, v in G.edges():
    #     res_original = G_copy[u][v].get('resistance')
    #     res_modified = G[u][v].get('resistance')
    #     if res_original != res_modified:
    #         differences.append(((u, v), res_original, res_modified))

    # if differences:
    #     print("Edges with changed resistance:")
    #     for edge, res_orig, res_mod in differences:
    #         print(f"Edge {edge}: {res_orig} -> {res_mod}")
    # else:
    #     print("No resistance changes detected.")
    # raise ValueError
    return G



