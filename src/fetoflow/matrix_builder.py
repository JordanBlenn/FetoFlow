import numpy as np
import scipy.sparse as sp
import networkx as nx


def create_matrices(G, bcs):
    """Create the full system matrices for the pressure-flow equations.
    
    Constructs the system matrices A and b that represent both the flow 
    conservation equations at nodes and the pressure-flow relationships on edges.
    The system handles both pressure and flow boundary conditions and multiple
    inlet nodes.

    Parameters
    ----------
    G : networkx.DiGraph
        Graph containing network topology and vessel properties
    bcs : dict
        Dictionary specifying boundary conditions with format:
        {
            "inlet": {"pressure": float or dict} or {"flow": float or dict},
            "outlet": {"pressure": float}
        }
        For multiple inlets, inlet values should be a dict mapping node IDs
        to boundary values.

    Returns
    -------
    scipy.sparse.csr_matrix
        System matrix A
    numpy.ndarray
        Right-hand side vector b

    Notes
    -----
    The system is constructed as:
    [A_flow     | A_elem] [p] = [b_flow]
    [A_pressure | A_res ] [q]   [b_pres]
    
    where:
    - A_flow: Flow conservation at nodes 
    - A_elem: Edge-node incidence
    - A_pressure: Node-edge incidence
    - A_res: Edge resistance terms
    - p: Node pressures
    - q: Edge flows
    - b_flow, b_pres: Boundary condition terms
    
    Matrix format is LIL during construction for efficient incremental
    assembly, then converted to CSR for computation.
    """
    # Get system dimensions
    n = G.number_of_nodes()  # Number of nodes in network
    m = G.number_of_edges()  # Number of edges (vessels)
    # Initialize system matrices in LIL format for efficient assembly
    A = sp.lil_matrix((n + m, n + m))  # System matrix: [n×n flow cons | n×m edge-node; m×n pressure | m×m resistance]
    b = np.zeros(n + m)                 # RHS vector: [n flow BCs; m pressure drops]

    # Determine boundary condition types from dict structure
    is_inlet_pressure = next(iter(bcs["inlet"].keys())) == "pressure"  # True if pressure BC at inlet
    is_inlet_flow = next(iter(bcs["inlet"].keys())) == "flow"         # True if flow BC at inlet
    is_outlet_pressure = next(iter(bcs["outlet"].keys())) == "pressure"  # Currently only pressure outlet supported
    # Validate boundary condition specification
    if not (is_inlet_pressure or is_inlet_flow):
        raise ValueError("No inlet boundary condition! Valid types are 'pressure' and 'flow'")
    if not is_outlet_pressure:
        raise ValueError("No outlet boundary condition! Valid types are 'pressure'")
    if is_inlet_pressure and is_inlet_flow:
        raise TypeError("Multiple inlet boundary conditions! Only 1 type (pressure or flow) allowed")

    # TODO MORE CHECKING VALID DICT TYPE???

    # Check for multiple inlet boundary conditions (different values at different inlets)
    multiple_inlet_bcs = False
    if is_inlet_pressure:
        inlet_bc = bcs["inlet"]["pressure"]
        if isinstance(inlet_bc, dict):  # Dict means different values per inlet
            multiple_inlet_bcs = True

    if is_inlet_flow:
        inlet_bc = bcs["inlet"]["flow"]
        if isinstance(inlet_bc, dict):  # Dict means different flows per inlet
            multiple_inlet_bcs = True

    # Get outlet pressure BC (currently only uniform pressure supported)
    outlet_bc = bcs["outlet"]["pressure"]  # TODO: Support different outlet pressures

    # Verify graph dimensions are consistent
    assert n == len(G.nodes), f"Number of nodes in the graph ({len(G.nodes)} does not match n ({n}).)"
    assert m == len(G.edges), f"Number of edges in the graph ({len(G.edges)} does not match m ({m}).)"

    # Assemble system starting with flow conservation equations
    for i in G.nodes:  # Loop over all nodes
        in_edges = G.in_edges(i)   # Edges entering node i
        out_edges = G.out_edges(i)  # Edges leaving node i
        
        # Handle inlet nodes (no incoming edges)
        if len(in_edges) == 0:  # TODO: Improve handling of multiple inlets/outlets
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


def create_small_matrices(G, bcs, branching_angles=False, non_linear_rheology=False):
    """Create reduced system matrices using boundary elimination.
    
    Constructs a reduced system by eliminating boundary nodes and forming
    a Schur complement system for the internal nodes. Optionally includes
    terms for bifurcation pressure losses and non-linear blood rheology.

    Parameters
    ----------
    G : networkx.DiGraph
        Graph containing network topology and vessel properties
    bcs : dict
        Boundary condition dictionary as in create_matrices()
    branching_angles : bool, optional
        Whether to include bifurcation pressure losses
    non_linear_rheology : bool, optional
        Whether to include non-linear viscosity effects

    Returns
    -------
    scipy.sparse.csr_matrix
        Reduced system matrix A 
    numpy.ndarray
        Reduced right-hand side vector b
    tuple
        Boundary condition information (bc_type, indices, values, inlet_idx)
    dict, optional
        Additional info for iterative solution if branching_angles=True:
        - branch_nodes: Dict mapping node indices to (entering, leaving) lists
        - branching_calc_matrices: Pre-computed matrices for efficiency 
        - branching_update_matrix: Initial resistance matrix
        - inlet_bc: BC type string
        - max_inlet_bc_flow: Maximum inlet flow value

    Notes
    -----
    The reduced system is formed by:
    1. Identifying boundary nodes from incidence matrix
    2. Forming Schur complement for internal nodes
    3. Incorporating nonlinear effects if requested
    4. Pre-computing matrices needed for iteration
    """
    # Initialize iteration options if needed
    if branching_angles:  # or any other nonlinear effects
        iter_options = {}
    else:
        iter_options = None
    # Validate boundary condition types
    is_inlet_pressure = next(iter(bcs["inlet"].keys())) == "pressure"  # Check inlet type
    is_inlet_flow = next(iter(bcs["inlet"].keys())) == "flow"
    is_outlet_pressure = next(iter(bcs["outlet"].keys())) == "pressure"  # Check outlet type
    
    # Ensure valid boundary conditions are specified
    if not (is_inlet_pressure or is_inlet_flow):
        raise ValueError("No inlet boundary condition! Valid types are 'pressure' and 'flow'")
    if not (is_outlet_pressure):
        raise ValueError("No outlet boundary condition! Valid types are 'pressure'")
    if is_inlet_pressure and is_inlet_flow:
        raise TypeError("Multiple inlet boundary conditions! Only 1 type allowed")
    else:
        bc_export = []  # Will store processed boundary condition info
        # Construct basic system matrices
        
        B = sp.csr_matrix(nx.incidence_matrix(G, oriented=True))  # CSR format for efficient row operations
        
        # Identify boundary nodes (degree 1) and internal nodes
        nnz_per_row = B.getnnz(axis=1)  # Count non-zero entries per row (node degree)
        boundary_indices = np.where(nnz_per_row == 1)[0]  # Nodes with one connection
        rest = np.where(~(nnz_per_row == 1))[0]  # Internal nodes
        
        # Initialize bifurcation tracking if needed
        if branching_angles:
            iter_options["branch_nodes"] = {}  
            # Loop through internal nodes looking for bifurcations
            for new_row, current_row in enumerate(rest):
                if B[current_row,:].nnz > 2:  # Node has >2 connections = bifurcation
                    # Get lists of entering and leaving vessel nodes
                    leaving_nodes = [n[1] for n in G.out_edges(current_row)] 
                    entering_nodes = [n[0] for n in G.in_edges(current_row)] 
                    
                    # Store bifurcation info: (entering, leaving, junction node)
                    iter_options["branch_nodes"][new_row] = (entering_nodes, leaving_nodes, current_row)
        else:
            iter_options = None  # No bifurcation tracking needed

        # Create W = diag(1/R)
        vals = np.array([1 / G[u][v]["resistance"] if G[u][v]["resistance"] != 0 else 0 
                        for u, v in G.edges()]) 
        W = sp.diags(vals, offsets=0).tocsc()  # Diagonal matrix of conductances
        
        # Pre-compute W*B^T for efficiency
        WBt = W @ B.transpose() 

        # Process boundary nodes and values
        boundary_vals = []  # Will store BC values in order
        outlet_idx = []    # Track outlet node indices
        inlet_idx = []     # Track inlet node indices
        error = False      # Flag for invalid flow BC specification
        
        # Loop through boundary nodes 
        for node in boundary_indices:
            # Check if inlet 
            if np.sum(B[node,:]) == -1:  # Inlet node
                if is_inlet_pressure:
                    if isinstance(bcs["inlet"]["pressure"], dict):
                        p_val =  bcs["inlet"]["pressure"].get(node)
                        if p_val:
                            boundary_vals.append(p_val)
                        else:
                            error = True
                    # Simple pressure BC - use same value
                    else:
                        boundary_vals.append(bcs["inlet"]["pressure"])
                else:
                    # Flow BC - handle both single and multiple values
                    if isinstance(bcs["inlet"]["flow"], dict):
                        # Multiple flow values - look up for this node
                        val = bcs["inlet"]["flow"].get(node)
                        if val:
                            boundary_vals.append(val)
                        else:
                            error = True  # Missing flow value for this inlet
                    else:
                        # Single flow value - use for all inlets
                        boundary_vals.append(bcs["inlet"]["flow"])
                inlet_idx.append(node)
            else:  # Outlet node
                outlet_idx.append(node)
                boundary_vals.append(bcs["outlet"]["pressure"])
        # Validate flow BC specification
        if error:
            raise ValueError(f"The values specified for the boundary inlets are not valid. "
                           f"The value should be the 'entering node' (the node you would use "
                           f"for a pressure inlet)\nThe valid nodes are: {np.array(inlet_idx)+1}")

        # Additional processing for flow boundary conditions
        if is_inlet_flow:
            edge_idx = []  # Track indices of inlet edges
            
            # Get list of inlet flow values
            if isinstance(bcs["inlet"]["flow"], dict):
                flow_in = list(bcs["inlet"]["flow"].values())  # Multiple specified flows
            else:
                # Single flow value replicated for each inlet
                flow_in = [bcs["inlet"]["flow"]]*(len(boundary_vals) - len(outlet_idx))
            
            # Find edges connected to inlet nodes
            for node in boundary_indices:
                if np.sum(B[node,:]) == -1:  # Inlet node
                    data = B[node,:].indices[0]  # Get connecting edge index
                    edge_idx.append(data)
                    
            # Create inverse resistance matrix for inlet edges
            a_inv_vals = 1 / vals[edge_idx]  # Get inverse resistances
            a_inv = sp.diags(a_inv_vals)      # Create diagonal matrix
            # Extract submatrices for weighted Laplacian
            Br = B[rest,:]                    # Reduced incidence matrix (internal nodes only)
            M = Br @ WBt[:,rest]              # Main diagonal block
            u = Br @ WBt[:,inlet_idx]         # Coupling to inlet nodes
            v = Br @ WBt[:,outlet_idx]        # Coupling to outlet nodes
            c = B[inlet_idx, :] @ WBt[:,rest] # Flow conservation at inlets
            
            # Create boundary value vectors
            p_out = sp.csc_array(np.array([bcs["outlet"]["pressure"]] * len(outlet_idx)).reshape(-1, 1))
            vp_out = v @ p_out                # Outlet pressure contribution
            qi = sp.csc_array(np.array(flow_in).reshape(-1, 1))  # Inlet flow vector
            
            # Form System
            A = M - u @ a_inv @ c             # System matrix with eliminated inlets
            b = -u @ a_inv @ qi - vp_out      # Modified RHS vector
            b = b.tocsc()                     # Convert to CSC for better multiplication

            # Package boundary condition info and matrices
            bc_export = ("Flow",              # BC type
                        np.array(boundary_indices),  # Boundary node indices 
                        np.array(boundary_vals),     # Boundary values
                        inlet_idx)                   # Inlet node indices
                        
            # Store matrices needed for bifurcation iteration
            branching_angles_matrices = [Br, u, c, qi, vp_out, edge_idx]
            
        else:
            # Simpler system for pressure boundary conditions
            Br = B[rest,:]  # Reduced incidence matrix
            A = Br @ WBt[:,rest]  # Direct Schur complement
            b = -Br @ WBt[:,boundary_indices] @ boundary_vals  # RHS with BCs
            
            # Package boundary info - no inlet distinction needed
            bc_export = ("Pressure",
                        np.array(boundary_indices),
                        np.array(boundary_vals),
                        None)  # No inlet tracking needed
                        
            # Only need reduced incidence for bifurcations
            branching_angles_matrices = [Br]

        # Package iteration options if needed
        if branching_angles:
            # Store matrices and info needed for bifurcation iteration
            iter_options["branching_calc_matrices"] = branching_angles_matrices  # Pre-computed matrices
            iter_options["branching_update_matrix"] = W  # Initial conductance matrix
            iter_options["inlet_bc"] = bc_export[0]     # BC type for scaling
            # Maximum inlet flow (0 for pressure BCs)
            iter_options["max_inlet_bc_flow"] = max(flow_in) if is_inlet_flow else 0
            return A, b, bc_export, iter_options
        
        # Basic return without iteration options
        return A, b, bc_export

