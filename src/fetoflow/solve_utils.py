import scipy.sparse as sps
import numpy as np
import time
from scipy.special import expit

from .matrix_builder import create_matrices
from .resistance_utils import calculate_resistance, calculate_viscosity_factor_from_radius
from .geometry_utils import update_geometry_with_pressures_and_flows

def solve_small_system(A, b, G,boundary_conditions,ill_conditioned=False,p0=None,current_p=None,max_iterations=None,restart=None):
    bc_type,boundary_indices,boundary_vals,inlet_idx = boundary_conditions

    p = sps.linalg.spsolve(A,b) # This solves for internal pressures!
    if ill_conditioned:
        internal_p = p.copy()
    q = np.zeros(shape=G.number_of_edges())
    if bc_type == "Flow":
        boundary_vals_current_iteration = boundary_vals.copy() # otherwise blowup!  
        for current_inlet in inlet_idx:
            adj_to_inlet = list(G.out_edges(current_inlet))[0][1]
            value_idx = np.where(boundary_indices == current_inlet)[0][0]
            p0 = boundary_vals[value_idx]*G[current_inlet][adj_to_inlet]["resistance"] + p[current_inlet]
            boundary_vals_current_iteration[value_idx] = p0
    else:
        boundary_vals_current_iteration = boundary_vals
    indices = np.argsort(boundary_indices)
    for idx,val in zip(boundary_indices[indices],boundary_vals_current_iteration[indices]):
        p = np.insert(p,idx,val)
    for u,v in G.edges():
        pu,pv = p[u],p[v]
        q[G[u][v]['edge_id']] = (pu-pv)/G[u][v]['resistance']
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    pressures = {node_id: p[node_id] for node_id in range(num_nodes)}
    flows = {edge_id: q[edge_id] for edge_id in range(num_edges)} #TODO: triple check these but they should be fine
    if ill_conditioned:
        return pressures,flows,internal_p
    return pressures, flows

def __solve_with_gmres(A, b, G,boundary_conditions,current_p=None):
    bc_type,boundary_indices,boundary_vals,inlet_idx = boundary_conditions
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

    
def update_small_matrix(G,diag_to_update,flows,iter_options):
    # depending on diverging/converging use associated edges to compute pressure drop coefficients for the current bifurcation/edge combo.
    # could look at only computing some of these once for efficiency, but i think since the flows need to be recalculated and updated anyway not worth
    t = time.time()
    if iter_options["branch_nodes"]:
        viscous_re_threshold = 100 # scale when Re < 100 -> otherwise the Mynard approx is not rly valid for low Re
        absolute_re_threshold = 300 # reynolds number should be at most 300 (based on umbilical artery/vein)
        rho = 1060 # blood density
        mu = 4e-3 # visc
        branching_dict = iter_options["branch_nodes"]
        for row in branching_dict.keys():
            # get required info
            entering_nodes,leaving_nodes,original_row = branching_dict[row]
            if len(leaving_nodes) > 1:
                datum_flow,datum_area = flows[G[entering_nodes[0]][original_row]["edge_id"]],np.pi*G[entering_nodes[0]][original_row]["radius"]**2
                max_flow = absolute_re_threshold*datum_area*mu/(rho*2*np.pi*G[entering_nodes[0]][original_row]["radius"])
                datum_flow = np.clip(datum_flow,1e-14,max_flow)
                for j in leaving_nodes:
                    flow_j,area_j = flows[G[original_row][j]["edge_id"]],np.pi*G[original_row][j]["radius"]**2
                    max_flow = absolute_re_threshold*area_j*mu/(rho*2*np.pi*G[original_row][j]["radius"])
                    flow_j = np.clip(flow_j,1e-14,max_flow)

                    coef = 1 - 1/((flow_j/datum_flow)*(datum_area/area_j))*np.cos((3/4*np.pi - 3/4*(np.pi - G[original_row][j]["bifurcation_angle"])))
                    p_loss = coef*rho*(flow_j/area_j)**2 + 1/2*rho*((datum_flow/datum_area)**2 - (flow_j/area_j)**2) # full pressure loss as per Mynard

                    Re = rho * (flow_j/area_j)* 2* G[original_row][j]["radius"]/mu
                    if Re < viscous_re_threshold:
                        p_loss*= Re/viscous_re_threshold

                    effective_res = p_loss/flow_j
                    diag_to_update[G[original_row][j]["edge_id"]] = 1/(1/diag_to_update[G[original_row][j]["edge_id"]] + effective_res)

            elif len(entering_nodes) > 1:
                datum_flow,datum_area = flows[G[original_row][leaving_nodes[0]]["edge_id"]],np.pi*G[original_row][leaving_nodes[0]]["radius"]**2
                max_flow = absolute_re_threshold*datum_area*mu/(rho*2*np.pi*G[entering_nodes[0]][original_row]["radius"])
                datum_flow = np.clip(datum_flow,1e-14,max_flow)
                p_loss = 0
                for j in entering_nodes:
                    flow_j,area_j = flows[G[j][original_row]["edge_id"]],np.pi*G[j][original_row]["radius"]**2
                    max_flow = absolute_re_threshold*area_j*mu/(rho*2*np.pi*G[j][original_row]["radius"])
                    flow_j = np.clip(flow_j,1e-14,max_flow)
                    coef = 1 - 1/((flow_j/datum_flow)*(datum_area/area_j))*np.cos((3/4*np.pi - 3/4*(np.pi - G[j][original_row]["bifurcation_angle"])))
                    p_loss += coef*rho*(flow_j/area_j)**2 + 1/2*rho*((datum_flow/datum_area)**2 - (flow_j/area_j)**2)
                if Re < viscous_re_threshold:
                    p_loss*= Re/viscous_re_threshold
                effective_res = p_loss/datum_flow
                diag_to_update[G[original_row][leaving_nodes[0]]["edge_id"]] = 1/(1/diag_to_update[G[original_row][leaving_nodes[0]]["edge_id"]] + effective_res)

    print(f"matrix update: {time.time() - t}")
    return diag_to_update




def solve_system(A, b, num_nodes, num_edges):
    x = sps.linalg.spsolve(A, b)
    pressures = {node_id: x[node_id] for node_id in range(num_nodes)}
    flows = {edge_id: x[num_nodes + edge_id] for edge_id in range(num_edges)}
    # TODO: Check this works
    return pressures, flows

def iterative_solve_small(A,b,G,bc_export,tol,info,alpha=1,maxiter=20,use_gmres=False,adaptive_stepping=True):
    s = time.time()
    W = info["branching_update_matrix"]
    if use_gmres:
        p,q,internal_p = solve_small_system(A,b,G,bc_export,ill_conditioned=True)
    else:
        p,q = solve_small_system(A,b,G,bc_export,ill_conditioned=False)
    p0 = np.array([p[node] for node in p.keys()])
    q0 = np.array([q[elem] for elem in q.keys()])
    print(f"Info |  Max Pressure: {np.max(p0)} | Min Pressure: {np.min(p0)} | Max flow: {np.max(q0)} | Min flow: {np.min(q0)}")

    p_mse = np.inf
    old_mse = np.inf
    iteration = 1
    print("Warm up period...")
    min_alpha = min(5*tol,1)
    max_alpha = 1
    flag = int(use_gmres)
    init_alpha = alpha
    if len(info["branching_calc_matrices"]) == 1:
        Br = info["branching_calc_matrices"][0]
    else:
        [Br,u_ainv_c] = info["branching_calc_matrices"]
    while p_mse > tol or flag:
        if iteration > maxiter:
            print("Maximum Iteration Count Reached. Consider constraining the learning rate parameter alpha (if observing cycling) or GMRES if MSE is very large")
            break
        diag_update = update_small_matrix(G,W.diagonal(),q,info)
        W_new = sps.diags(diag_update,offsets=0).tocsc()
        if len(info["branching_calc_matrices"]) == 1:
            A_new = Br @ W_new @ Br.T
        else:
            A_new = Br @ W_new @ Br.T - u_ainv_c
        if use_gmres:
            p,q,internal_p,flag = __solve_with_gmres(A_new,b,G,bc_export,current_p=internal_p)
        else:
            p,q = solve_small_system(A_new,b,G,bc_export,ill_conditioned=False)
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


def solve_iterative_system(G, A, b, num_nodes, num_edges, bcs, viscosity_model, mu, capillary_model, capillary_parameters, tol=0.01, max_solve_time=120): # max solve time in seconds
    #TODO: Make solve iterative small system work as well
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

def update_A_matrix(A, pressures, flows, G,n,m, flow_dependent_viscosity=False, branching_angles=False, elastic_vessels=False,):
    #TODO: Write this function
    if branching_angles:
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


def iterative_solve(A,b,G,bc_export,c0=0):
    p,q = solve_small_system(A,b,G,bc_export)
