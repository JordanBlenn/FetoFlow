import numpy as np


def calculate_resistance(
    G,
    viscosity_model="constant",
    mu=0.33600e-02,
    capillary_model="analytical2015",
    capillary_parameters=None,
):
    if capillary_parameters is None:
        capillary_parameters = {
            "num_series": 3,
            "num_parallel": 6,
            "num_generations": 3,
            "num_convolutes_per_branch": 10,
            "capillary_hematocrit": 0.45,
            "total_capillary_length": 3.0 / 1000,
            "capillary_convolute_radius": 7.2e-6,
            "radius_cap_artery": 1.5e-5,
            "radius_cap_scaling_factor": 2,
            "pries_mu": 0.4e-2,
            "total_segment_length": 1.5 / 1000,
            "capillary_length": 3.0 / 6 / 1000,
            "segment_length": 1.5 / 10 / 1000,
        }

    valid_viscocity_models = {"constant", "pries_network", "pries_vessel", "flow_dependent"} # Flow-dependent is also pries_vessel but with flow nonlinear stuff.

    if viscosity_model not in valid_viscocity_models:
        raise ValueError(
            f"Invalid Viscocity Model: {viscosity_model}. \
                         Valid viscosity models are 'constant' (default), 'pries_network', 'pries_vessel', and 'flow_dependent'."
        )

    viscosity_factor = 1 # Should be for non-capillary vessels always
    for u, v in G.edges():
        # Specific mu for blood vessel with flow_dependent viscosity
        if viscosity_model == "flow_dependent":
            mu = G[u][v]['mu']
            viscosity_factor = G[u][v]['viscosity_factor']
        
        if not G[u][v]["vessel_type"] == "capillary_equivalent":
            # Non-linear blood rheology doesn't apply to any non-capillary vessels as they are too large.
            # Flow-dependent it should - but viscosity factor should be 1 if not changed.
            G[u][v]["resistance"] = 8 * mu * viscosity_factor * G[u][v]["length"] / (np.pi * G[u][v]["radius"] ** 4)
        else:  # capillary resistance
            radius_in_artery = list(G.in_edges(u, data=True))[0][-1][
                "radius"
            ]  # assuming 1 vessel into the terminal arterial node, use that. We check this earlier.
            radius_in_vein = list(G.out_edges(v, data=True))[0][-1][
                "radius"
            ]  # assuming 1 vessel out of the terminal venous node, use that. We check this earlier.
            G[u][v]["resistance"] = calculate_capillary_equivalent_resistance(
                radius_in_artery=radius_in_artery,
                radius_in_vein=radius_in_vein,
                viscosity_model=viscosity_model,
                mu=mu,
                capillary_model=capillary_model,
                capillary_parameters=capillary_parameters,
            )

    return G


def calculate_viscosity_factor_from_radius(radius, hematocrit=0.45):
    # TODO pries_vessel vs pries_network? Or is this not important here and only in capillary stuff?
    # THIS IS DIRECTLY COPIED FROM REPROSIM CODE SO I HOPE IT IS RIGHT.
    radius_um = radius * 1_000_000  # This formula's constants are based on micrometers I believe.TODO CHECK THE PAPER
    # print("Radius um: " + str(radius_um))
    control_hematocrit = 0.45
    beta = 4.0 / (1.0 + np.exp(-0.0593 * (2.0 * radius_um - 6.74)))
    viscosity_factor = (
        1.0
        + (np.exp(hematocrit * beta) - 1.0)
        / (np.exp(control_hematocrit * beta) - 1.0)
        * (110 * np.exp(-2.848 * radius_um) + 3.0 - 3.45 * np.exp(-0.07 * radius_um))
    ) / 4.0
    return viscosity_factor


def calculate_convolute_resistance(
    n_series=3, n_parallel=6, mu=0.33600e-02, visc_factor=1, capillary_length=0.0005, capillary_convolute_radius=7.2e-6
):
    # TODO: Check this is valid. FROM REPROSIM CODE.
    R_cap = (
        (8.0 * mu * visc_factor * capillary_length) / (np.pi * capillary_convolute_radius**4.0) * n_series / n_parallel
    )  # !%resistance of each capillary convolute segment
    return R_cap


def calculate_capillary_equivalent_resistance(
    radius_in_artery,
    radius_in_vein,
    viscosity_model="constant",
    mu=0.33600e-02,
    capillary_model="analytical2015",
    capillary_parameters=None,
):
    if capillary_parameters is None:
        capillary_parameters = {
            "num_series": 3,
            "num_parallel": 6,
            "num_generations": 3,
            "num_convolutes_per_branch": 10,
            "capillary_hematocrit": 0.45,
            "total_capillary_length": 3.0 / 1000,
            "capillary_convolute_radius": 7.2e-6,
            "radius_cap_artery": 1.5e-5,
            "radius_cap_scaling_factor": 2,
            "pries_mu": 0.4e-2,
            "total_segment_length": 1.5 / 1000,
            "capillary_length": 3.0 / 6 / 1000, #TODO: Check this should be divided by num_parallel, I think it should be.
            "segment_length": 1.5 / 10 / 1000,
        }

    if viscosity_model == "constant":
        viscosity_factor = 1
        convolute_resistance = calculate_convolute_resistance(
            mu=mu,
            capillary_length=capillary_parameters["capillary_length"],
            capillary_convolute_radius=capillary_parameters["capillary_convolute_radius"],
        )  # mu also is default here
    elif viscosity_model == "pries_vessel":
        mu = capillary_parameters["pries_mu"]
        # TODO: Documentation between cap convolute radius (each convolute) and rad_cap_artery - final intermediate villous segment!
        convolute_resistance = calculate_convolute_resistance(
            mu=mu,
            visc_factor=calculate_viscosity_factor_from_radius(
                radius=capillary_parameters["capillary_convolute_radius"], hematocrit=capillary_parameters["capillary_hematocrit"]
            ),
            capillary_length=capillary_parameters["capillary_length"],
            capillary_convolute_radius=capillary_parameters["capillary_convolute_radius"],
        )
    elif viscosity_model == "pries_network":
        mu = capillary_parameters["pries_mu"]
        viscosity_factor = 1
        convolute_resistance = calculate_convolute_resistance(
            mu=mu,
            capillary_length=capillary_parameters["capillary_length"],
            capillary_convolute_radius=capillary_parameters["capillary_convolute_radius"],
        )
    elif viscosity_model == "flow_dependent": # Fuck me I need to do this per individual element?
        mu = capillary_parameters["pries_mu"]
        viscosity_factor = 1
        convolute_resistance = calculate_convolute_resistance(
            mu=mu,
            capillary_length=capillary_parameters["capillary_length"],
            capillary_convolute_radius=capillary_parameters["capillary_convolute_radius"],
        )
    else:
        raise NotImplementedError("Other Viscosities not implemented!")

    # Calculate resistance of a single convolute
    R_c = convolute_resistance
    # print("Convolute Resistance: " + str(convolute_resistance))

    # Vein capillary radius is scaled from the artery capillary radius
    radius_cap_artery = capillary_parameters["radius_cap_artery"]
    radius_cap_vein = radius_cap_artery * capillary_parameters["radius_cap_scaling_factor"]
    num_generations = capillary_parameters["num_generations"]  # TODO: Tidy up, do I want all variables like this, or keep as indexing dict?
    # Radii along the tree of intermediate villi are linearly interpolated from the capillary radius (smallest, included as radius of nth generation)
    # to the radius of the incoming artery/vein (largest, not included in the interpolation)
    artery_radii = np.linspace(radius_in_artery, radius_cap_artery, num_generations, endpoint=False)
    vein_radii = np.linspace(radius_in_vein, radius_cap_vein, num_generations, endpoint=False)
    # Reverse order of radii for bottom up calculations
    artery_radii = artery_radii[::-1]
    vein_radii = vein_radii[::-1]

    valid_model_types = ["analytical2015"]

    # Calculate effective resistance of each group of convolutes.
    # TODO: WE WILL NEED GOOD DOCUMENTATION EXPLAINING HOW THIS WORKS.
    if capillary_model == "analytical2015":
        # Loop through generations bottom-up:
        for height, (artery_radius, vein_radius) in enumerate(zip(artery_radii, vein_radii)):  # height = num_gerations - current_gen
            # Calculate segment resistances
            if viscosity_model == "pries_vessel":
                viscosity_factor = calculate_viscosity_factor_from_radius(artery_radius)

            segment_resistance_artery = 8 * mu * viscosity_factor * capillary_parameters["segment_length"] / (np.pi * artery_radius**4)

            if viscosity_model == "pries_vessel":
                viscosity_factor = calculate_viscosity_factor_from_radius(vein_radius)

            segment_resistance_vein = 8 * mu * viscosity_factor * capillary_parameters["segment_length"] / (np.pi * vein_radius**4)

            for __ in range(capillary_parameters["num_convolutes_per_branch"] - 1):
                # Add in series with segment resistances of veins and arteries, then in parallel with another convolute
                R_s = segment_resistance_artery + segment_resistance_vein + R_c
                R_c = 1 / (1 / R_s + 1 / convolute_resistance)
            R_level = segment_resistance_artery + segment_resistance_vein + R_c
            # Final iteration (first generation) should not add another convolute in parallel as this is where the tree initially splits
            R_c = 1 / (2 / R_level + 1 / convolute_resistance * (not height == num_generations - 1))
    else:
        raise ValueError(f"Invalid capillary model type: {capillary_model}.\n Current valid model types: {valid_model_types}")
    return R_c