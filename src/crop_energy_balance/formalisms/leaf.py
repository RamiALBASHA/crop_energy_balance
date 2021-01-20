def calc_leaf_boundary_conductance(wind_speed_at_canopy_height: float,
                                   characteristic_length: float = 0.01,
                                   shape_parameter: float = 0.01) -> float:
    """Calculates bulk boundary layer conductance (for both sides of leaves) at the scale of an individual leaf.

    Args:
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter

    Returns:
        [m h-1] bulk boundary layer conductance (for both sides of leaves) at the scale of an individual leaf

    """
    return 3600 * shape_parameter * (wind_speed_at_canopy_height / characteristic_length) ** 0.5


def calc_stomatal_sensibility(air_vapor_pressure_deficit: float,
                              shape_parameter: float) -> float:
    """Calculates the effect of air vapor pressure deficit on stomatal aperture.

    Args:
        air_vapor_pressure_deficit: [kPa] vapor pressure deficit of the air at canopy source height
        shape_parameter: [kPa] empirical shape parameter

    Returns:
        [-] stomatal closure fraction due to elevated air vapor pressure deficit

    """
    assert (shape_parameter != 0), 'The value of `shape_parameter` must be greater than zero.'
    return 1.0 / (1.0 + air_vapor_pressure_deficit / shape_parameter)


def calc_stomatal_conductance(residual_stomatal_conductance: float,
                              maximum_stomatal_conductance: float,
                              absorbed_irradiance: float,
                              shape_parameter: float,
                              stomatal_sensibility_to_water_status: float) -> float:
    """Calculates stomatal conductivity to water vapor for a unit leaf surface area.

    Args:
        residual_stomatal_conductance: [m h-1] residual (minimum) stomatal conductance
        maximum_stomatal_conductance: [m h-1] maximum stomatal conductance
        absorbed_irradiance: [W m-2leaf] absorbed photosynthetically active radiation (PAR)
        shape_parameter: [W m-2leaf] an empirical parameter that regulates the response of stomatal conductance to the
            absorbed PAR
        stomatal_sensibility_to_water_status: [-] stomatal closure fraction due to water stress

    Returns:
        [m h-1] stomatal conductivity to water vapor for a unit leaf surface area
    """
    stomatal_sensibility_to_absorbed_irradiance = (absorbed_irradiance / (shape_parameter + absorbed_irradiance))
    return residual_stomatal_conductance + maximum_stomatal_conductance * (
            stomatal_sensibility_to_absorbed_irradiance * stomatal_sensibility_to_water_status)
