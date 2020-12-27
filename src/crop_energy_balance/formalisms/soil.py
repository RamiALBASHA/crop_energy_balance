from math import exp, log


def calc_eddy_diffusivity(wind_speed: float,
                          canopy_height: float,
                          measurement_height: float,
                          zero_displacement_height: float,
                          canopy_roughness_length_for_momentum: float,
                          von_karman_constant: float) -> float:
    """Calculates the value of eddy diffusivity of water vapor and heat transfer at canopy height.

    Args:
        wind_speed: [m h-1] wind speed at reference height
        canopy_height: [m] average height of the canopy
        measurement_height: [m] height at which wind speed in measured
        zero_displacement_height: [m] zero displacement height
        canopy_roughness_length_for_momentum: [m] roughness length for momentum transfer of the canopy
        von_karman_constant: [-] von Karman constant

    Returns:
        [m2 h-1] eddy diffusivity at canopy height

    References:
        Choudhury and Monteith, 1988
            A four-layer model for the heat budget of homogeneous land surfaces.
            Q. J. Roy. Meteor. Soc. 114, 373 – 398.
        Montes et al. 2014
            A three-source SVAT modeling of evaporation: Application to the seasonal dynamics of a grassed vineyard.
            Agricultural and Forest Meteorology 191, 64 – 80.
            Eq. C5
    """
    return (von_karman_constant ** 2 * wind_speed * (canopy_height - zero_displacement_height)) / (
        log((measurement_height - zero_displacement_height) / canopy_roughness_length_for_momentum))


def calc_boundary_resistance(wind_speed: float,
                             canopy_height: float,
                             measurement_height: float,
                             zero_displacement_height: float,
                             canopy_roughness_length_for_momentum: float,
                             soil_roughness_length_for_momentum: float = 0.01,
                             shape_parameter: float = 2.5,
                             von_karman_constant: float = 0.41) -> float:
    """Calculates the bulk aerodynamic resistance resistance between the soil surface and the source height.

    Args:
        wind_speed: [m h-1] wind speed at a given hour
        canopy_height: [m] average canopy height
        measurement_height: [m] height at which meteorological measurements are made
        zero_displacement_height: [m] zero displacement height
        canopy_roughness_length_for_momentum: [m] roughness length for momentum transfer of the canopy
        soil_roughness_length_for_momentum: [m] roughness length for momentum transfer of bare soil
        shape_parameter: [-] a shape parameter
        von_karman_constant: [-] von Karman constant

    Returns:
        [h m-1] bulk soil boundary layer resistance

    References:
        Choudhury and Monteith, 1988
            A four-layer model for the heat budget of homogeneous land surfaces.
            Q. J. Roy. Meteor. Soc. 114, 373 – 398.
        Montes et al. 2014
            A three-source SVAT modeling of evaporation: Application to the seasonal dynamics of a grassed vineyard.
            Agricultural and Forest Meteorology 191, 64 – 80.
            Eq. C4
    """

    eddy_diffusivity = calc_eddy_diffusivity(
        wind_speed=wind_speed,
        canopy_height=canopy_height,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        canopy_roughness_length_for_momentum=canopy_roughness_length_for_momentum,
        von_karman_constant=von_karman_constant)

    scaling_factor = exp(-shape_parameter * soil_roughness_length_for_momentum / canopy_height) - exp(
        -shape_parameter * (zero_displacement_height + canopy_roughness_length_for_momentum) / canopy_height)

    return canopy_height * exp(shape_parameter) / (shape_parameter * eddy_diffusivity) * scaling_factor


def calc_surface_resistance(soil_saturation_ratio: float,
                            shape_parameter_1: float = 8.206,
                            shape_parameter_2: float = 4.255) -> float:
    """Calculates the bulk soil surface resistance.

    Args:
        soil_saturation_ratio: [-] ratio of actual to potential volumetric water content in the soil
        shape_parameter_1: [s m-1] empirical shape parameter
        shape_parameter_2: [s m-1] empirical shape parameter

    Returns:
        [h m-1] bulk soil surface resistance to water vapor transfer

    References:
        Sellers et al., 1992
            Relations between surface conductance and spectral vegetation indices at intermediaite (100 m2 to 15 km2)
            length scales.
            Journal of Geophysical Research 97, 19,033 - 19,059
    """

    return 1.0 / 3600.0 * exp(shape_parameter_1 - shape_parameter_2 * soil_saturation_ratio)


def calc_heat_flux(net_above_ground_radiation: float,
                   is_diurnal: bool) -> float:
    """Calculates the net heat flux density into the soil substrates.

    Args:
        net_above_ground_radiation: [W m-2ground] the net above-ground radiation
        is_diurnal: `True` during the daylight hours, otherwise `False`

    Returns:
        [W m-2ground] the net heat flux density into the soil substrates

    """
    return 0.1 * net_above_ground_radiation if is_diurnal else 0.5 * net_above_ground_radiation


def calc_net_longwave_radiation(canopy_top_net_longwave_radiation: float,
                                canopy_leaf_area_index: float,
                                diffuse_black_extinction_coefficient: float) -> float:
    """Calculates net long wave radiation exchange at the soil surface.

    Args:
        canopy_top_net_longwave_radiation: [W m-2ground] net long wave radiation at the top of the canopy
        canopy_leaf_area_index: [m2leaf m-2ground] total leaf area index of the canopy
        diffuse_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of diffuse irradiance
            through a canopy of black leaves

    Returns:
        [W m-2ground]: net long wave radiation exchange at the soil surface

    References:
        Leuning et al. 1995
            Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies.
            Plant, Cell and Environment 18, 1183 - 1200.
    """
    scaling_factor = exp(-diffuse_black_extinction_coefficient * canopy_leaf_area_index)
    return canopy_top_net_longwave_radiation * scaling_factor
