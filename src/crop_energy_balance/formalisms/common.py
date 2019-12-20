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


def calc_soil_heat_flux(net_above_ground_radiation: float,
                        is_diurnal: bool) -> float:
    """Calculates the net heat flux density into the soil substrates.

    Args:
        net_above_ground_radiation: [W m-2ground] the net above-ground radiation
        is_diurnal: `True` during the daylight hours, otherwise `False`

    Returns:
        [W m-2ground] the net heat flux density into the soil substrates

    """
    return 0.1 * net_above_ground_radiation if is_diurnal else 0.5 * net_above_ground_radiation


def calc_component_net_radiation(net_shortwave_radiation: float,
                                 net_longwave_radiation: float,
                                 is_soil: bool) -> float:
    """Calculates the net radiation flux density of a canopy component per unit ground surface area.

    Args:
        net_shortwave_radiation: [W m-2ground] net shortwave flux density radiation of the component
        net_longwave_radiation: [W m-2ground] net longwave flux density radiation of the component
        is_soil: `True` if the calculation is performed for the soil component, otherwise `False`

    Returns:
        [W m-2ground] net radiation flux density of the canopy component per unit ground surface area

    """
    if not is_soil:
        return net_shortwave_radiation + net_longwave_radiation
    else:
        net_above_ground_radiation = net_shortwave_radiation + net_longwave_radiation
        return net_above_ground_radiation - calc_soil_heat_flux(net_above_ground_radiation, net_shortwave_radiation > 0)


def calc_canopy_aerodynamic_resistance(wind_speed: float,
                                       canopy_height: float,
                                       reference_height: float,
                                       von_karman_cst: float) -> float:
    """Calculates air resistance to water vapor transfer between the source height and the reference height.

    Args:
        wind_speed: [m h-1] wind speed at measurement height
        canopy_height: [m] average height of the canopy
        reference_height: [m] height at which wind speed in measured
        von_karman_cst: [-] von Karman constant

    Returns:
        [h m-1] air resistance to water vapor transfer between the source height and the reference height
    """
    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)
    zero_displacement_height_ = calc_zero_displacement_height(canopy_height)
    canopy_roughness_length_for_momentum_ = calc_canopy_roughness_length_for_momentum(canopy_height)
    canopy_roughness_length_for_heat_transfer_ = calc_canopy_roughness_length_for_heat_transfer(canopy_height)

    return 1.0 / (wind_speed * von_karman_cst ** 2) * (
            log((reference_height - zero_displacement_height_) / canopy_roughness_length_for_momentum_) *
            log((reference_height - zero_displacement_height_) / canopy_roughness_length_for_heat_transfer_))


def calc_vapor_pressure_slope(temperature: float) -> float:
    """Calculates the slope of air vapor pressure curve at a given air temperature.

    Args:
        temperature: [°C] air temperature

    Returns:
        [kPa °K-1] the slope of vapor pressure curve at the given air temperature

    """
    return 4098 * (0.6108 * exp((17.27 * temperature) / (temperature + 237.3))) / ((temperature + 237.3) ** 2)


def calc_penman_evaporative_energy(canopy_aerodynamic_resistance: float,
                                   canopy_net_radiation: float,
                                   vapor_pressure_slope: float,
                                   vapor_pressure_deficit: float,
                                   psychrometric_constant: float,
                                   air_density: float,
                                   air_specific_heat_capacity: float) -> float:
    """Calculates the evapotranspiration energy flux density according to Penman's formula.

    Args:
        canopy_aerodynamic_resistance: [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        canopy_net_radiation: [W m-2ground] net radiation flux density of the entire canopy per unit ground surface area
        vapor_pressure_slope: [kPa K-1] the slope of vapor pressure curve at a given air temperature
        vapor_pressure_deficit: [kPa] vapor pressure deficit at the reference height
        psychrometric_constant: [kPa K-1] the psychrometric constant
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [W m-2ground] the evapotranspiration energy flux density according to Penman's formula

    """
    # canopy_net_radiation = max(0.0, canopy_net_radiation)  # to avoid calculating negative transpiration fluxes
    return (vapor_pressure_slope * canopy_net_radiation + (
            air_density * air_specific_heat_capacity * vapor_pressure_deficit) / canopy_aerodynamic_resistance) / (
                   vapor_pressure_slope + psychrometric_constant)


def calc_canopy_lumped_aerodynamic_resistance(canopy_aerodynamic_resistance: float,
                                              vapor_pressure_slope: float,
                                              psychrometric_constant: float) -> float:
    """Calculates Lhomme's lumped aerodynamic resistance (R0).

    Args:
        canopy_aerodynamic_resistance: [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        vapor_pressure_slope: [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_constant: [kPa K-1] the psychrometric constant

    Returns:
        [h m-1] Lhomme's lumped aerodynamic resistance (R0)

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 7
    """
    return (1 + vapor_pressure_slope / psychrometric_constant) * canopy_aerodynamic_resistance


def calc_component_composed_resistance(surface_resistance: float,
                                       boundary_layer_resistance: float,
                                       vapor_pressure_slope: float,
                                       psychrometric_constant: float,
                                       stomatal_density_factor: int) -> float:
    """Calculates Lhomme's lumped boundary and surface resistance (Ri) for a given canopy component.

    Args:
        surface_resistance: [h m-1] surface resistance of the component
        boundary_layer_resistance: [h m-1] boundary layer resistance of the component
        vapor_pressure_slope: [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_constant: [kPa K-1] the psychrometric constant
        stomatal_density_factor: [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    Returns:
        [h m-1] Lhomme's lumped boundary and surface resistance (Ri) for the given canopy component

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 8
    """
    return surface_resistance + boundary_layer_resistance * (
            stomatal_density_factor + vapor_pressure_slope / psychrometric_constant)


def calc_component_composed_conductance(composed_boundary_and_surface_resistance: float,
                                        sum_composed_boundary_and_surface_conductances: float,
                                        canopy_lumped_aerodynamic_resistance: float) -> float:
    """Calculates Lhomme's composed boundary and surface conductance (Pi) of a given canopy component

    Args:
        composed_boundary_and_surface_resistance: [h m-1] Lhomme's lumped boundary and surface resistance for
            the given component
        sum_composed_boundary_and_surface_conductances: [m h-1] the sum of Lhomme's lumped boundary and surface
            resistance over all canopy components
        canopy_lumped_aerodynamic_resistance: [h m-1] Lhomme's lumped aerodynamic resistance

    Returns:
        [m h-1] Lhomme's composed boundary and surface conductance (Pi) of the given canopy component

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 13
    """
    return 1.0 / (composed_boundary_and_surface_resistance * (
            1.0 + canopy_lumped_aerodynamic_resistance * sum_composed_boundary_and_surface_conductances))


def calc_penman_monteith_evaporative_energy(components_indices: list,
                                            canopy_lumped_aerodynamic_resistance: float,
                                            penman_evaporative_energy: float,
                                            composed_boundary_and_surface_conductances: dict,
                                            net_radiation_fluxes: dict,
                                            boundary_layer_resistances: dict,
                                            vapor_pressure_slope: float,
                                            psychrometric_constant: float) -> float:
    """Calculates canopy evaporative energy.

    Args:
        components_indices: [-] list of dictionary keys of all components of the canopy
        canopy_lumped_aerodynamic_resistance: [h m-1] Lhomme's lumped aerodynamic resistance (R0)
        penman_evaporative_energy: [W m-2ground] the evapotranspiration energy flux density according to
            Penman's formula
        composed_boundary_and_surface_conductances: [m h-1] Lhomme's composed boundary and surface conductance
            (Pi) of all canopy components
        net_radiation_fluxes: [W m-2ground] net radiation flux densities of all canopy components
        boundary_layer_resistances: [h m-1] boundary layer reistance to heat transfer of all canopy components
        vapor_pressure_slope: [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_constant: [kPa K-1] the psychrometric constant

    Returns:
        [W m-2ground] canopy evaporative energy

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 12
    """
    sum_p = 0.0
    sum_p_a_r = 0.0
    for i in components_indices:
        sum_p += composed_boundary_and_surface_conductances[i]
        sum_p_a_r += (
                composed_boundary_and_surface_conductances[i] * net_radiation_fluxes[i] * boundary_layer_resistances[i])
    return canopy_lumped_aerodynamic_resistance * penman_evaporative_energy * sum_p + (
            vapor_pressure_slope / psychrometric_constant) * sum_p_a_r


def calc_component_evaporative_energy(net_radiation: float,
                                      boundary_layer_resistance: float,
                                      lumped_boundary_and_surface_resistance: float,
                                      canopy_lumped_aerodynamic_resistance: float,
                                      penman_evaporative_energy: float,
                                      penman_monteith_evaporative_energy: float,
                                      vapor_pressure_slope: float,
                                      psychrometric_constant: float) -> float:
    """Calculates the evaporative energy for a given canopy component.

    Args:
        net_radiation: [W m-2ground] net radiation flux density of the given component
        boundary_layer_resistance: [h m-1] component boundary layer resistance to heat transfer
        lumped_boundary_and_surface_resistance: [h m-1] Lhomme's lumped boundary and surface resistance (Ri)
            for the given component
        canopy_lumped_aerodynamic_resistance: [h m-1] Lhomme's lumped aerodynamic resistance (R0)
        penman_evaporative_energy: [W m-2ground] the evapotranspiration energy flux density according to
            Penman's formula
        penman_monteith_evaporative_energy: [W m-2ground] canopy evaporative energy
        vapor_pressure_slope: [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_constant: [kPa K-1] the psychrometric constant

    Returns:
        [W m-2ground] the evaporative energy for the given canopy component

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 11
    """
    energy_driven_evaporation = canopy_lumped_aerodynamic_resistance * (
            penman_evaporative_energy - penman_monteith_evaporative_energy)
    radiation_driven_evaporation = boundary_layer_resistance * net_radiation * (
            vapor_pressure_slope / psychrometric_constant)
    return 1.0 / lumped_boundary_and_surface_resistance * (energy_driven_evaporation + radiation_driven_evaporation)


def calc_canopy_temperature(air_temperature: float,
                            canopy_aerodynamic_resistance: float,
                            canopy_net_radiation: float,
                            penman_monteith_evaporative_energy: float,
                            air_density: float,
                            air_specific_heat_capacity: float) -> float:
    """Calculates air temperature at source height.

    Args:
        air_temperature: [K] air temperature at measurement height
        canopy_aerodynamic_resistance: [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        canopy_net_radiation: [W m-2ground] canopy net radiation
        penman_monteith_evaporative_energy: [W m-2ground] canopy evaporative energy
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [K] air temperature at source height
    """
    return air_temperature + (canopy_aerodynamic_resistance / (air_density * air_specific_heat_capacity)) * (
            canopy_net_radiation - penman_monteith_evaporative_energy)


def calc_component_temperature(canopy_temperature: float,
                               boundary_layer_resistance: float,
                               component_net_radiation: float,
                               component_evaporative_energy: float,
                               air_density: float,
                               air_specific_heat_capacity: float) -> float:
    """Calculates air temperature at source height.

    Args:
        canopy_temperature: [K] air temperature at source height
        boundary_layer_resistance: [h m-1] resistance to heat transfer at the given component
        component_net_radiation: [W m-2ground] canopy net radiation
        component_evaporative_energy: [W m-2ground] canopy evaporative energy
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [K] air temperature at source height
    """
    return canopy_temperature + (boundary_layer_resistance / (air_density * air_specific_heat_capacity)) * (
            component_net_radiation - component_evaporative_energy)


def calc_temperature_step(previous_value: float,
                          actual_value: float,
                          step_fraction: float = 0.5):
    """Calculates the temperature step value between two consecutive energy balance calculations.

    Args:
        previous_value: [K] previous calculated temperature
        actual_value: [K] actual calculated temperature
        step_fraction: [-] fraction of the entire step (`actual_value - previous_value`) to be used

    Returns:
        [K] the temperature step value between two consecutive energy balance calculations
    """
    return step_fraction * (actual_value - previous_value)


def calc_component_net_longwave_radiation(air_temperature: float,
                                          lower_leaf_area_index: float,
                                          atmoshperic_emissivity: float,
                                          extinction_coef: float,
                                          stefan_boltzman_constant: float) -> float:
    """Calculates net long wave radiation of a canopy component.

    Args:
        air_temperature: [K] air temperature
        lower_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        atmoshperic_emissivity: [-] atmospheric emissivity to longwave radiation
        extinction_coef: [m2ground m-2leaf] extinction coefficient of longwave radiation, assumed equal to
            the extinction coefficient of the diffuse photosynthetically active radiation for black leaves
        stefan_boltzman_constant: [W m-2 K-4] Stefan-Boltzmann constant

    Returns:
        net long wave radiation of a canopy component
    """
    scaling_factor = 1.0 / extinction_coef * exp(-extinction_coef * lower_leaf_area_index)
    return - (1 - atmoshperic_emissivity) * stefan_boltzman_constant * air_temperature ** 4 * scaling_factor
