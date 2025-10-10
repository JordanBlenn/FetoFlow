import networkx as nx
import numpy as np


def create_geometry(
    nodes,
    elements,
    inlet_radius,
    strahler_ratio_arteries,
    arteries_only=False,
    outlet_vein_radius=None,
    strahler_ratio_veins=None,
    fields=None,
    default_mu=0.33600e-02,
    default_hematocrit=0.45,
):
    G = nx.DiGraph()
    num_terminal_arterial_nodes = 0
    # Read in Nodes and Edges
    for node_id, coordinates in nodes.items():
        G.add_node(node_id, x=coordinates[0], y=coordinates[1], z=coordinates[2])
    for edge_id, (node_from, node_to) in enumerate(elements):
        if fields:
            res = fields.get("resistance")
            if res:
                res = fields.get(edge_id, 0)
            else:
                res = 0
            # .get returns None by default if not found
            radius = fields.get("radius")
            if radius:
                radius = radius.get(edge_id)
            else:
                radius = None
        else:
            res = 0.0
            radius = None
        length = calcLength(G, node_from, node_to)
        G.add_edge(
            node_from,
            node_to,
            edge_id=edge_id,
            resistance=res,
            length=length,
            radius=radius,
            strahler=None,
            vessel_type="artery",
            bifurcation_angle=0,
            mu=default_mu,
            hematocrit=default_hematocrit,
            viscosity_factor=1,
        )

    # Find all input nodes (to ensure we give every element a strahler ordering)
    input_nodes = []
    for node in G.nodes():
        if len(G.in_edges(node)) == 0:
            input_nodes.append(node)
            # print(node)
    # Add artery radii via strahler ordering
    max_strahler = 0
    for input_node in input_nodes:
        for u, v in G.out_edges(input_node):
            G = update_strahlers(G, u, v)  # update for each input node should work
            max_strahler = max(max_strahler, G[u][v]["strahler"])

    if fields and fields.get("radius"):
        for u, v in G.edges():
            if G[u][v]["radius"] is None:
                elem_strahler = G[u][v]["strahler"]
                # need to update R according to the last specified radius in the subtree
                radius_found = False
                inlet_radius_updated = inlet_radius
                out_node = u
                sub_tree_strahler = max_strahler
                while not radius_found:
                    edges_to_check = list(G.in_edges(out_node))
                    # should only be one edge coming in
                    if len(edges_to_check) == 0:
                        radius_found = True  # no previously set radii in vessel's predecessors
                    else:
                        in_node = edges_to_check[0][0]
                        current_rad = G[in_node][out_node]["radius"]
                        if current_rad is not None and current_rad != 0:
                            inlet_radius_updated = current_rad
                            radius_found = True
                            sub_tree_strahler = G[in_node][out_node]["strahler"]
                        else:
                            out_node = in_node
                G[u][v]["radius"] = inlet_radius_updated * strahler_ratio_arteries ** (elem_strahler - sub_tree_strahler)
    else:
        for u, v in G.edges():
            elem_strahler = G[u][v]["strahler"]
            G[u][v]["radius"] = inlet_radius * strahler_ratio_arteries ** (elem_strahler - max_strahler)

    # Create the venous mesh
    if not arteries_only:
        num_artery_nodes = G.number_of_nodes()  # use for scaling to keep numeric based indexing
        num_artery_edges = G.number_of_edges()
        # Get terminal nodes
        terminal_nodes = [node for node, out_degree in G.out_degree() if out_degree == 0]
        num_terminal_arterial_nodes = len(terminal_nodes)
        venous_mesh = create_venous_mesh(
            G,
            num_artery_nodes,
            num_artery_edges,
            num_terminal_arterial_nodes,
            outlet_vein_radius,
            strahler_ratio_veins,
            max_strahler,
        )
        # list of artery terminal node indices
        # add venous mesh to graph
        assert max(G.nodes) < min(venous_mesh.nodes), "Venous mesh node ids overlap with arterial node ids."
        G = nx.compose(G, venous_mesh)
        # Add edges to each terminal node with equivalent capillary network resistance
        edge_id_tracker = len(elements)  # Use this to easily increment edge_id correctly for the new edges I am adding.
        for i, arterial_node in enumerate(terminal_nodes):
            venous_node = arterial_node + num_artery_nodes

            assert len(G.in_edges(arterial_node)) == 1, f"Terminal artery node has {len(G.in_edges(arterial_node))} entering it, should be 1 only."
            assert len(G.out_edges(venous_node)) == 1, f"Terminal vein node has {len(G.out_edges(venous_node))} exiting it, should be 1 only."

            G.add_edge(
                arterial_node,
                venous_node,
                edge_id=edge_id_tracker + i,
                resistance=None,
                length=None,
                radius=0.0,
                strahler=0.0,
                vessel_type="capillary_equivalent",
                mu=default_mu,
                hematocrit=default_hematocrit,
                viscosity_factor=1,

            )
            # Update length of these vessels - we calculate the resistance in 'calculate_resistance()' function.
            G[arterial_node][venous_node]["length"] = calcLength(G, arterial_node, venous_node)
            # Leave radius as 0 - this allows visualisations not including the capillary networks,
            #  which would only show a single vessel anyway (not the tree of the intermediate and terminal villi).

    return G


def create_venous_mesh(
    G,
    num_artery_nodes,
    num_artery_edges,
    num_terminal_arterial_nodes,
    outlet_vein_radius,
    strahler_ratio_veins,
    max_strahler,
):
    venous_mesh = G.copy()
    # first n/2 nodes = artery nodes going down.
    # Next n/2 nodes = vein nodes in same order artery nodes were (i.e. likely inlet first, getting smaller as we go).
    nx.relabel_nodes(venous_mesh, lambda node_id: node_id + num_artery_nodes, copy=False)  # 0-based indexing works here
    # Update radii
    for u, v in venous_mesh.edges():
        # If anastomosis, remove edge as it doesn't exist in the veins
        if venous_mesh[u][v]["vessel_type"] == "anastomosis":
            G.remove_edge(u, v)
            continue

        # Edge ids for veins are after capillaries
        venous_mesh[u][v]["edge_id"] = venous_mesh[u][v]["edge_id"] + num_artery_edges + num_terminal_arterial_nodes
        venous_mesh[u][v]["vessel_type"] = "vein"
        # Strahler ordering for veins
        elemStrahler = venous_mesh[u][v]["strahler"]
        venous_mesh[u][v]["radius"] = outlet_vein_radius * strahler_ratio_veins ** (elemStrahler - max_strahler)
    # Reverse the direction of all arcs
    venous_mesh = venous_mesh.reverse()

    return venous_mesh


def update_strahlers(G, node_in, node_out):
    # Input the input node(s) -  or call on each input node. Updates the strahler field.
    # This is a recursive function and will be slow for now.
    child_edges = G.out_edges(node_out)
    # # Base case: No children, strahler == 1.
    if len(child_edges) == 0:
        G[node_in][node_out]["strahler"] = 1
        return G

    max_child_strahler = 0
    max_child_strahler_count = 0

    for u, v in child_edges:
        if G[u][v]["strahler"] is None:
            G = update_strahlers(G, u, v)  # Returns graph object with strahler updated.

        if G[u][v]["strahler"] > max_child_strahler:
            max_child_strahler = G[u][v]["strahler"]
            max_child_strahler_count = 1
        elif G[u][v]["strahler"] == max_child_strahler:
            max_child_strahler_count += 1

    assert max_child_strahler_count > 0
    # Actually assign this arc's strahler now.
    if max_child_strahler_count == 1:
        G[node_in][node_out]["strahler"] = max_child_strahler
    else:
        G[node_in][node_out]["strahler"] = max_child_strahler + 1  # 2 arcs coming in with same max value

    return G


def calcLength(G, u, v):
    return np.sqrt(np.sum([(G.nodes[u][coord] - G.nodes[v][coord]) ** 2 for coord in ["x", "y", "z"]])) / 1000  # mm to m!
    # TODO: Fix unit conversions and makr work for anything etc. i.e. specify units somewhere as an input argument at the start

def create_anastomosis(G, node_from, node_to, radius=None, mu=0.33600e-02):
    # NOTE HERE: RADIUS IS IN mm!!!!!
    # Todo: make sure this is clear.
    # TODO: DIGRAPH STUFF???. No we can just have a negative flow along the anastomosis.
    # TODO: probably write this as 2 separate functions, this one which is in here with the graph stuff, and one which is user-facing and calls this one
    u = node_from - 1
    v = node_to - 1  # Update from 1- to 0-based indexing\

    if u not in G:
        raise ValueError(
            f"Node {node_from} (ipnode indexing)/Node {u} (networkX indexing) does not exist in the networkX graph. Perhaps you need to call create_geometry() first?"
        )
    if v not in G:
        raise ValueError(
            f"Node {node_to} (ipnode indexing)/Node {v} (networkX indexing) does not exist in the networkX graph. Perhaps you need to call create_geometry() first?"
        )
    if u == v:
        raise ValueError(f"Anastomosis cannot connect the same node to itself. Node number is {node_from} (ipnode indexing)/{u} (networkX indexing).")

    # only exists as aterial connection
    G.add_edge(
        u,
        v,
        edge_id=G.number_of_edges(),
        resistance=0.0,
        length=None,
        radius=0.0,
        strahler=0.0,
        vessel_type="anastomosis",
        mu=mu,
        hematocrit=0.45,# TODO PARAMETERISE
        viscosity_factor=1
    )
    # Old implementation:
    # - defines a radius of the anastomosis which is used
    # Our alternative:
    # - Provide a warning if this happens, but just use the maxiumum radii of the vessels leaving the nodes connecting the anastomosis.
    # TODO: Confirm this is fine implementation

    G[u][v]["length"] = calcLength(G, node_from, node_to)

    # For strahler, just take the max of any child strahler.
    max_child_radius = 0.0
    max_child_strahler = 0.0
    for node in (u, v):
        for child_u, child_v in G.out_edges(node):
            child_radius = G[child_u][child_v]["radius"]
            child_strahler = G[child_u][child_v]["strahler"]
            if not child_radius or not child_strahler:
                raise ValueError("Strahler values and radii have not been set yet. make sure you do this before creating an anastomosis.")
            max_child_radius = max(max_child_radius, child_radius)
            max_child_strahler = max(max_child_strahler, child_strahler)

    if radius:
        if not (isinstance(radius, int) or isinstance(radius, float)):
            raise ValueError(f"Hyrtl anastomosis radius is invalid type {type(radius)}. Valid types are float or int")

        G[u][v]["radius"] = radius / 1000  # mm to m!
    else:
        G[u][v]["radius"] = max_child_radius
    # Strahler
    G[u][v]["strahler"] = max_child_strahler  # The code will previously break if strahlers have not been already set.
    # Resistance calculation: #TODO make sure it updates properly if we have other viscosities etc.
    # Note: if calcualte_resistance() is called after this function, the result will be overwritten. This is probably the order we want:
    # - calculate_geometry()
    # - create_venous_mesh()
    # - create_anatsomosis()
    # - calculate_resistance()
    G[u][v]["resistance"] = 8 * mu * G[u][v]["length"] / (np.pi * G[u][v]["radius"] ** 4)

    return G


def update_geometry_with_pressures_and_flows(G, pressures, flows):
    for node_id in G.nodes():
        G.nodes[node_id]["pressure"] = pressures.loc[node_id]["pressure"]
    for u, v in G.edges():
        G[u][v]["flow"] = flows.loc[G[u][v]["edge_id"]]["flow"]
    return G

def calculate_branching_angles(G):
    # double check...
    for n in G.nodes():
        out_edges = list(G.out_edges(n))
        in_edges = list(G.in_edges(n))
        if len(out_edges) > 1:
            for __,out_node in out_edges:
                out_vec = np.array([G.nodes[n][coord] - G.nodes[out_node][coord] for coord in ["x","y","z"]])
                out_norm= out_vec/np.linalg.norm(out_vec)
                in_vec = np.array([G.nodes[in_edges[0][0]][coord] - G.nodes[n][coord] for coord in ["x","y","z"]])
                in_norm = in_vec/np.linalg.norm(in_vec)
                dot = np.clip(in_norm @ out_norm,-1,1)
                G[n][out_node]["bifurcation_angle"] = np.pi - np.arccos(dot) #phi_j from Mynard
                
        elif len(in_edges) > 1:
            for in_node,__ in in_edges:
                in_vec = np.array([G.nodes[in_node][coord] - G.nodes[n][coord]  for coord in ["x","y","z"]])
                in_norm = in_vec/np.linalg.norm(in_vec)
                out_vec = np.array([G.nodes[n][coord] - G.nodes[out_edges[0][1]][coord] for coord in ["x","y","z"]])
                out_norm = out_vec/np.linalg.norm(out_vec)
                dot = np.clip(in_norm @ out_norm,-1,1)
                G[in_node][n]["bifurcation_angle"] = np.pi - np.arccos(dot) #phi_j from Mynard
        else:
            continue
    return

        


