import numpy as np
import scipy.sparse as sp
import networkx as nx


def create_matrices(G, n, m, bcs):  # n = num nodes, m = num elements. NOT SIZE OF ARRAY AS IT DYNAMICALLY RESIZES NOW.
    # Start as LiL matrix
    A = sp.lil_matrix((n + m, n + m))
    b = np.zeros(n + m)

    # Check for Flow or pressure BC's
    is_inlet_pressure = next(iter(bcs["inlet"].keys())) == "pressure"
    is_inlet_flow = next(iter(bcs["inlet"].keys())) == "flow"
    is_outlet_pressure = next(iter(bcs["outlet"].keys())) == "pressure"
    if not (is_inlet_pressure or is_inlet_flow):
        raise ValueError("No inlet boundary condition! Valid types are 'pressure' and 'flow'")
    if not (is_outlet_pressure):
        raise ValueError("No outlet boundary condition! Valid types are 'pressure'")  # Make parameterised?
    if is_inlet_pressure and is_inlet_flow:
        raise TypeError("Multiple  inlet boundary conditions! Only 1 type allowed")

    # TODO MORE CHECKING VALID DICT TYPE???

    # Check for single or multiple inlet BCs
    multiple_inlet_bcs = False
    if is_inlet_pressure:
        inlet_bc = bcs["inlet"]["pressure"]
        if isinstance(inlet_bc, dict):
            multiple_inlet_bcs = True

    if is_inlet_flow:
        inlet_bc = bcs["inlet"]["flow"]
        if isinstance(inlet_bc, dict):
            multiple_inlet_bcs = True

    # TODO Other outlet bcs?
    outlet_bc = bcs["outlet"]["pressure"]

    assert n == len(G.nodes), f"Number of nodes in the graph ({len(G.nodes)} does not match n ({n}).)"
    assert m == len(G.edges), f"Number of edges in the graph ({len(G.edges)} does not match m ({m}).)"

    # Add Flow Equations first
    for i in G.nodes:  # i is the numeric node id
        in_edges = G.in_edges(i)
        out_edges = G.out_edges(i)
        if (
            len(in_edges) == 0
        ):  # Input node TODO Better handling around multiple input and output nodes?? This only looks at the graph object currently, doesn't check the dict.
            if is_inlet_pressure:
                A[i, i] = 1
                if not multiple_inlet_bcs:
                    b[i] = inlet_bc
                else:
                    if i not in inlet_bc.keys():
                        raise KeyError(
                            f"Node {i + 1} is not defined as an inlet node when specifying pressure boundary conditions."
                            f"Note: this refers to Node {i + 1} from the ipnode file but Node {i} in the NetworkX graph due to 1- vs 0-based indexing."
                        )
                    b[i] = inlet_bc[i]
            else:  # Flow BC
                n_rows = A.shape[0]
                n_cols = A.shape[1]
                # Create a dummy node connecting to the input node with no pressure difference between them (i.e. a resistance of 0).
                # Define the flow between them.
                dummy_rows = sp.lil_matrix((2, n_cols + 2))  # Dummy rows get appended at bottom of matrix: [dummy_node; dummy_arc]
                dummy_rows[0, -2] = 1  # Pressure at dummy node
                dummy_rows[0, i] = -1  # Pressure at inlet node
                dummy_rows[1, -1] = 1  # Dummy arc value - known inlet flow

                A.resize((n_rows + 2, n_cols + 2))  # Add extra rows to A Matrix.
                A[n_rows, :] = dummy_rows[0, :]
                A[n_rows + 1, :] = dummy_rows[1, :]

                if not multiple_inlet_bcs:
                    b = np.append(b, [0, inlet_bc])
                else:
                    if i not in inlet_bc.keys():
                        raise KeyError(
                            f"Element {i + 1} is not defined as an inlet node when specifying flow boundary conditions."
                            f"Note: this refers to Element {i + 1} from the ipelem file but Element {i} in the NetworkX graph due to 1- vs 0-based indexing."
                        )
                    b = np.append(b, [0, inlet_bc[i]])
                # Flow in = flow out for inlet node.
                A[i, -1] = 1  # Inlet flow is known
                for u, v in out_edges:
                    index = G[u][v]["edge_id"]
                    A[i, n + index] = -1  # Refers to existing element(s), not dummy node

        elif len(out_edges) == 0:  # outlet node
            A[i, i] = 1
            b[i] = outlet_bc  # Assumes all outlet BC's are the same for terminal nodes. Outlet Pressure

        # If not terminal, setup matrix for equation
        else:
            # Not a terminal node
            for u, v in in_edges:
                index = G[u][v]["edge_id"]
                A[i, n + index] = 1  # CHECK THIS
            for u, v in out_edges:
                index = G[u][v]["edge_id"]
                A[i, n + index] = -1

    # Add pressure Equations
    for u, v in G.edges():  # U and V are the node from id and node_to id.
        edge_id = G[u][v]["edge_id"]
        R = G[u][v]["resistance"]
        A[edge_id + n, u] = 1
        A[edge_id + n, v] = -1
        A[edge_id + n, edge_id + n] = -R  # This fixes the previous bug of using i and enumerate
    # Convert to CSR
    A = A.tocsr()

    return A, b


def create_small_matrices(G, bcs,branching_angles=False,non_linear_rheology=False):
    if branching_angles: # or any other kwarg:
        iter_options = {}
    else:
        iter_options = None
    # Check for Flow or pressure BC's
    is_inlet_pressure = next(iter(bcs["inlet"].keys())) == "pressure"
    is_inlet_flow = next(iter(bcs["inlet"].keys())) == "flow"
    is_outlet_pressure = next(iter(bcs["outlet"].keys())) == "pressure"
    if not (is_inlet_pressure or is_inlet_flow):
        raise ValueError("No inlet boundary condition! Valid types are 'pressure' and 'flow'")
    if not (is_outlet_pressure):
        raise ValueError("No outlet boundary condition! Valid types are 'pressure'")
    if is_inlet_pressure and is_inlet_flow:
        raise TypeError("Multiple  inlet boundary conditions! Only 1 type allowed")
    else:
        bc_export = []
        # Construct required matrices
        B = sp.csr_matrix(nx.incidence_matrix(G,oriented=True)) # needs to be csr for indexing inlet/outlet effectively
        nnz_per_row = B.getnnz(axis=1)
        boundary_indices = np.where(nnz_per_row == 1)[0]
        rest = np.where(~(nnz_per_row == 1))[0]
        if branching_angles:
            iter_options["branch_nodes"] = {}
            for new_row,current_row in enumerate(rest):
                if B[current_row,:].nnz > 2: # otherwise no branching angle effect
                    leaving_nodes = [n[1] for n in G.out_edges(current_row)]  
                    entering_nodes = [n[0] for n in G.in_edges(current_row)]
                    iter_options["branch_nodes"][new_row] = (entering_nodes,leaving_nodes,current_row) # new row will be the associated index in the r x r matrix A
        else:
            iter_options = None

        vals = np.array([1 / G[u][v]["resistance"] if G[u][v]["resistance"] != 0 else 0 for u, v in G.edges()])
        W = sp.diags(vals,offsets=0).tocsc()
        WBt = W @ B.transpose() # precompute for efficiency

        ## TODO lot of this done in bc_utils. Will clean up once implementation complete.
        boundary_vals = []
        outlet_idx = []
        inlet_idx = []
        for node in boundary_indices:
            if np.sum(B[node,:]) == -1:
                boundary_vals.append(bcs["inlet"]["pressure"] if is_inlet_pressure else bcs["inlet"]["flow"])
                inlet_idx.append(node)
            else:
                outlet_idx.append(node)
                boundary_vals.append(bcs["outlet"]["pressure"])

        if is_inlet_flow:
            edge_idx = []
            if isinstance(bcs["inlet"]["flow"],dict):
                edge_idx = bcs["inlet"]["flow"].keys()
                flow_in = bcs["inlet"]["flow"].values()
            else:
                flow_in = [bcs["inlet"]["flow"]]*(len(boundary_vals) - len(outlet_idx))
                for node in boundary_indices:
                    if np.sum(B[node,:]) == -1:
                        data = B[node,:].indices[0]
                        edge_idx.append(data)
            vals = W.diagonal()
            a_inv_vals = 1 / vals[edge_idx]
            a_inv = sp.diags(a_inv_vals)
            Br = B[rest, :]
            M = Br @ WBt[:,rest]
            u = Br @  WBt[:,inlet_idx]
            v = Br @ WBt[:,outlet_idx]
            c = B[inlet_idx, :] @  WBt[:,rest]

            u_ainv_c = u @ a_inv @ c
                
            qi = sp.csc_array(np.array(flow_in).reshape(-1, 1))
            p_out = sp.csc_array(np.array([bcs["outlet"]["pressure"]] * len(outlet_idx)).reshape(-1, 1))
            A = M - u_ainv_c
            b = -u @ a_inv @ qi - v @ p_out
            b = b.tocsc() # more efficient for sparse multiplication

            bc_export = ("Flow",np.array(boundary_indices),np.array(boundary_vals),inlet_idx)
            branching_angles_matrices = [Br,u_ainv_c]
        else:
            # Pressure BCs
            Br = B[rest,:]
            A = Br @ WBt[:,rest]
            b = -Br @ WBt[:,boundary_indices] @ boundary_vals
            bc_export = ("Pressure",np.array(boundary_indices),np.array(boundary_vals),None)
            branching_angles_matrices = [Br]

        # A = A.tocsc()
        if branching_angles:
            iter_options["branching_calc_matrices"] = branching_angles_matrices
            iter_options["branching_update_matrix"] = W
        return A, b, bc_export, iter_options

