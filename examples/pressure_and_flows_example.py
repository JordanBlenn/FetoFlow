from fetoflow.pressure_flow_utils import pressures_and_flows
import pandas as pd
def main():

    node_filename = "sample_geometry/3_branchlen_tree.ipnode"
    element_filename = "sample_geometry/3_branchlen_tree.ipelem"
    boundary_conditions = {
        "inlet_pressure" : 6650,
        "outlet_pressure" : 2660
    }
    inlet_radius = 1.8/1000
    strahler_ratio_arteries = 1.38
    
    # call pressure and flows function
    G = pressures_and_flows(
        node_filename,
        element_filename,
        boundary_conditions,
        inlet_radius,
        strahler_ratio_arteries,
        input_directory=".",
        output_directory="./output_data",
        flow_output_filename="flow_values.csv",
        pressure_output_filename="pressure_values.csv",
        arteries_only=False,
        viscosity_model="constant",
        vessel_type="rigid",
        outlet_vein_radius=4.0,
        strahler_ratio_veins=1.46,
        anastomosis=None,
        mu=0.33600e-2,  # This is the non-capillary viscosity value used
        capillary_model="analytical2015",
        capillary_parameters=None, 
        radius_filename=None,
        other_field_filenames=None,  
        verbose=False,
        time_statistics=False,
        return_graph=False,
    )

if __name__ == '__main__':
    main()
