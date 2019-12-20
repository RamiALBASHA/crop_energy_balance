from math import exp

from crop_energy_balance.formalisms import common


def calc_leaf_layer_boundary_conductance(wind_speed_at_canopy_height: float,
                                         upper_leaf_area_index: float,
                                         lower_leaf_area_index: float,
                                         wind_speed_extinction_coefficient: float = 0.5,
                                         characteristic_length: float = 0.01,
                                         shape_parameter: float = 0.01) -> float:
    """Calculates bulk layer boundary layer conductance (for both sides of leaves).

    Args:
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        upper_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        wind_speed_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter

    Returns:
        [m h-1] Calculates bulk layer boundary layer conductance (for both sides of leaves)

    """
    leaf_boundary_conductance = common.calc_leaf_boundary_conductance(wind_speed_at_canopy_height,
                                                                      characteristic_length,
                                                                      shape_parameter)
    scaling_factor = 2.0 / wind_speed_extinction_coefficient * (
            exp(-0.5 * wind_speed_extinction_coefficient * upper_leaf_area_index) -
            exp(-0.5 * wind_speed_extinction_coefficient * lower_leaf_area_index))

    return leaf_boundary_conductance * scaling_factor


def calc_leaf_layer_boundary_resistance_to_vapor(wind_speed_at_canopy_height: float,
                                                 upper_leaf_area_index: float,
                                                 lower_leaf_area_index: float,
                                                 wind_speed_extinction_coefficient: float,
                                                 characteristic_length: float,
                                                 shape_parameter: float,
                                                 stomatal_density_factor: int) -> float:
    """Calculates the bulk leaf layer resistance to water vapor transfer.

    Args:
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        upper_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        wind_speed_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-1/2] an empirical shape parameter
        stomatal_density_factor: [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    Returns:
        [h m-1] bulk leaf layer resistance to water vapor transfer

    """
    boundary_conductance = calc_leaf_layer_boundary_conductance(wind_speed_at_canopy_height,
                                                                upper_leaf_area_index,
                                                                lower_leaf_area_index,
                                                                wind_speed_extinction_coefficient,
                                                                characteristic_length,
                                                                shape_parameter)

    return 1.0 / (boundary_conductance / stomatal_density_factor)


def calc_leaf_layer_surface_conductance_to_vapor(absorbed_irradiance: float,
                                                 upper_cumulative_leaf_area_index: float,
                                                 lower_cumulative_leaf_area_index: float,
                                                 stomatal_sensibility_to_water_status: float,
                                                 global_extinction_coefficient: float,
                                                 maximum_stomatal_conductance: float,
                                                 residual_stomatal_conductance: float,
                                                 shape_parameter: float = 105) -> float:
    """Calculates the bulk surface conductance of a leaf layer.

    Args:
        absorbed_irradiance: [W m-2ground] absorbed photosynthetically active radiation
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index above the considered layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index below the considered layer
        stomatal_sensibility_to_water_status: [-] stomatal closure fraction due to water stress
        global_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of the entire PAR wavelength into
            the canopy
        maximum_stomatal_conductance: [m h-1] maximum stomatal conductance
        residual_stomatal_conductance: [m h-1] residual (minimum) stomatal conductance
        shape_parameter: [W m-2leaf] an empirical parameter to regulate the shape of stomatal conductance response
            to absorbed PAR

    Returns:
        [m h-1] bulk surface conductance of the leaf layer
    """
    scaling_factor = 1.0 / global_extinction_coefficient * log(
        (global_extinction_coefficient * absorbed_irradiance *
         exp(-global_extinction_coefficient * upper_cumulative_leaf_area_index) + shape_parameter) /
        (global_extinction_coefficient * absorbed_irradiance *
         exp(-global_extinction_coefficient * lower_cumulative_leaf_area_index) + shape_parameter))

    return residual_stomatal_conductance * (
            1 - scaling_factor) + maximum_stomatal_conductance * stomatal_sensibility_to_water_status * scaling_factor
