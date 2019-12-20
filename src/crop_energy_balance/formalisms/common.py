from math import log, exp, sqrt


def calc_psychrometric_constant(atmospheric_pressure: float,
                                air_specific_heat_capacity: float = 2.8e-4,
                                latent_heat_for_vaporization: float = 0.678,
                                vapor_to_dry_air_molecular_weight: float = 0.622) -> float:
    """Calculates the psychrometric constant.

    Args:
        atmospheric_pressure: [kPa] atmospheric pressure
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure
        latent_heat_for_vaporization: [w h g-1] latent heat for vaporization
        vapor_to_dry_air_molecular_weight: [-] ratio of the molecular weights of water vapor to dry air

    Returns:
        [kPa K-1] the psychrometric constant under the given atmospheric pressure

    References:
        Allen et al. 1998
            FAO Irrigation and Drainage Paper No. 56.
            Eq. 8
    """
    return air_specific_heat_capacity * atmospheric_pressure / (
            vapor_to_dry_air_molecular_weight * latent_heat_for_vaporization)


def calc_atmospheric_emissivity(air_vapor_pressure: float,
                                air_temperature: float) -> float:
    """Calculates the atmospheric emissivity to thermal infrared waves.

    Args:
        air_vapor_pressure: [kPa] air vapor pressure
        air_temperature: [K] air temperature

    Returns:
        [-] atmospheric emissivity
    """

    return 1.24 * (10. * air_vapor_pressure / air_temperature) ** (1. / 7.)


def calc_zero_displacement_height(canopy_height: float) -> float:
    """Calculates zero displacement height of water vapor between the crop and the atmosphere

    Args:
        canopy_height: [m] average height of the canopy

    Returns:
        [m] zero displacement height of water vapor between the crop and the atmosphere
    """
    return 0.67 * canopy_height


def calc_canopy_roughness_length_for_momentum(canopy_height: float) -> float:
    """Calculates roughness length for momentum

    Args:
        canopy_height: [m] average height of the canopy

    Returns:
        [m] roughness length for momentum
    """
    return 0.123 * canopy_height


def calc_canopy_roughness_length_for_heat_transfer(canopy_height: float) -> float:
    """Calculates roughness length for heat transfer

    Args:
        canopy_height: [m] average height of the canopy

    Returns:
        [m] roughness length for heat transfer
    """
    return 0.0123 * canopy_height


def calc_wind_speed_at_canopy_height(wind_speed: float,
                                     canopy_height: float,
                                     measurement_height: float) -> float:
    """Calculates hourly wind speed values at canopy height

    Args:
        wind_speed: [m h-1] hourly wind speed at measurement height
        canopy_height: [m] height of the canopy
        measurement_height: [m] height at which meteorological measurements are made

    Returns:
        [m h-1] wind speed at canopy height
    """

    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)
    d = calc_zero_displacement_height(canopy_height)
    z0u = calc_canopy_roughness_length_for_momentum(canopy_height)

    return max(10.0e-10, wind_speed * log((canopy_height - d) / z0u) / log((measurement_height - d) / z0u))


def calc_turbulent_diffusivity(von_karman_constant: float,
                               wind_speed: float,
                               canopy_height: float,
                               zero_displacement_height: float,
                               canopy_roughness_length_for_momentum: float,
                               measurement_height: float) -> float:
    """Calculates the turbulent (eddy) diffusivity of water vapor at canopy height

    Args:
        von_karman_constant: [-] von Karman constant
        wind_speed: [m h-1] wind speed at reference height
        canopy_height: [m] average height of the canopy
        zero_displacement_height: [m] zero displacement height
        canopy_roughness_length_for_momentum: [m] roughness length for momentum
        measurement_height: [m] height at which wind speed in measured

    Returns:
        [m2 h-1] turbulent (eddy) diffusivity of water vapor at canopy height
    """
    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)

    return (von_karman_constant ** 2 * wind_speed * (canopy_height - zero_displacement_height)) / (
        log((measurement_height - zero_displacement_height) / canopy_roughness_length_for_momentum))


def calc_soil_boundary_resistance(canopy_height: float,
                                  wind_speed: float,
                                  measurement_height: float,
                                  soil_roughness_length_for_momentum: float = 0.01,
                                  shape_parameter: float = 2.5,
                                  von_karman_constant: float = 0.41) -> float:
    """Calculates the bulk soil boundary layer resistance.

    Args:
        canopy_height: [m] average canopy height
        wind_speed: [m h-1] wind speed at a given hour
        measurement_height: [m] height at which meteorological measurements are made
        soil_roughness_length_for_momentum: [m] soil roughness length for momentum
        shape_parameter: [-] a shape parameter
        von_karman_constant: [-] von Karman constant

    Returns:
        [h m-1] bulk soil boundary layer resistance

    References:
        Choudhury and Monteith, 1988
            A four-layer model for the heat budget of homogeneous land surfaces.
            Q. J. Roy. Meteor. Soc. 114, 373 – 398.
    """

    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)
    zero_plane_displacement_height = calc_zero_displacement_height(canopy_height)
    canopy_roughness_length_for_momentum = calc_canopy_roughness_length_for_momentum(canopy_height)

    eddy_diffusivity = calc_turbulent_diffusivity(von_karman_constant,
                                                  wind_speed,
                                                  canopy_height,
                                                  zero_plane_displacement_height,
                                                  canopy_roughness_length_for_momentum,
                                                  measurement_height)

    scaling_factor = \
        exp(-shape_parameter * soil_roughness_length_for_momentum / canopy_height) - \
        exp(-shape_parameter * (
                zero_plane_displacement_height + canopy_roughness_length_for_momentum) / canopy_height)

    return canopy_height * exp(shape_parameter) / (shape_parameter * eddy_diffusivity) * scaling_factor


def calc_soil_surface_resistance(soil_saturation_ratio: float,
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
    return 3600 * shape_parameter * sqrt(wind_speed_at_canopy_height / characteristic_length)


def calc_leaf_layer_boundary_resistance_to_heat(leaf_layer_boundary_conductance: float) -> float:
    """Calculates the bulk leaf layer resistance to heat transfer.

    Args:
        leaf_layer_boundary_conductance: [m h-1] Calculates bulk layer boundary layer conductance (for both blade sides)

    Returns:
        [h m-1] bulk leaf layer resistance to heat transfer

    """
    return 1.0 / max(1.e-6, leaf_layer_boundary_conductance)


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

