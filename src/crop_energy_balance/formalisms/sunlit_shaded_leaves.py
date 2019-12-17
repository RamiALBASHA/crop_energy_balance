from math import exp, log, sqrt

from crop_energy_balance.formalisms.common import (calc_zero_displacement_height,
                                                   calc_canopy_roughness_length_for_momentum,
                                                   calc_canopy_roughness_length_for_heat_transfer)


def calc_leaf_layer_boundary_conductance(wind_speed_at_canopy_height, upper_leaf_area_index, lower_leaf_area_index,
                                         wind_speed_extinction_coef=0.5, characteristic_length=0.01,
                                         shape_parameter=0.01):
    """Calculates bulk layer boundary layer conductance (for both sides of leaves).

    Args:
        wind_speed_at_canopy_height (float): [m s-1] local wind speed in the vicinity of the leaf
        upper_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        wind_speed_extinction_coef (float): [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length (float): [m] characteristic leaf length in the direction of the wind
        shape_parameter (float): [m s-0.5] an empirical shape parameter

    Returns:
        (float): [m h-1] Calculates bulk layer boundary layer conductance (for both sides of leaves)

    """
    return 3600 * 2.0 * shape_parameter / wind_speed_extinction_coef * sqrt(
        wind_speed_at_canopy_height / characteristic_length) * (
                   exp(-0.5 * wind_speed_extinction_coef * upper_leaf_area_index) -
                   exp(-0.5 * wind_speed_extinction_coef * lower_leaf_area_index))


def calc_leaf_layer_boundary_resistance_to_vapor(wind_speed_at_canopy_height, upper_leaf_area_index,
                                                 lower_leaf_area_index,
                                                 wind_speed_extinction_coef, characteristic_length, shape_parameter,
                                                 stomatal_density_factor):
    """Calculates the bulk leaf layer resistance to water vapor transfer.

    Args:
        wind_speed_at_canopy_height (float): [m s-1] local wind speed in the vicinity of the leaf
        upper_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        wind_speed_extinction_coef (float): [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length (float): [m] characteristic leaf length in the direction of the wind
        shape_parameter (float): [m s-1/2] an empirical shape parameter
        stomatal_density_factor (int): [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    Returns:
        (float): [h m-1] bulk leaf layer resistance to water vapor transfer

    """
    boundary_conductance = calc_leaf_layer_boundary_conductance(wind_speed_at_canopy_height,
                                                                upper_leaf_area_index,
                                                                lower_leaf_area_index,
                                                                wind_speed_extinction_coef,
                                                                characteristic_length,
                                                                shape_parameter)

    return 1.0 / (boundary_conductance / stomatal_density_factor)


def calc_leaf_layer_boundary_layer_resistance_to_heat(leaf_layer_boundary_layer_conductance):
    """Calculates the bulk leaf layer resistance to heat transfer.

    Args:
        leaf_layer_boundary_layer_conductance (float): [m h-1] Calculates bulk layer boundary layer conductance
            (for both sides of leaves)

    Returns:
        (float): [h m-1] bulk leaf layer resistance to heat transfer

    """
    return 1.0 / leaf_layer_boundary_layer_conductance


def calc_stomatal_sensibility(air_vapor_pressure_deficit, shape_parameter):
    """Calculates the effect of air vapor pressure deficit on stomatal aperture.

    Args:
        air_vapor_pressure_deficit (float): [kPa] vapor pressure deficit of the air at canopy source height
        shape_parameter (float): [kPa] empirical shape parameter

    Returns:
        (float): [-] stomatal closure fraction due to elevated air vapor pressure deficit

    """
    assert (shape_parameter != 0), 'The value of `shape_parameter` must be greater than zero.'
    return 1.0 / (1.0 + air_vapor_pressure_deficit / shape_parameter)


def calc_leaf_layer_surface_conductance_to_vapor(incident_par, maximum_stomatal_conductance, stomatal_sensibility,
                                                 upper_cumulative_leaf_area_index, lower_cumulative_leaf_area_index,
                                                 absorbed_par_50=105, global_extinction_coefficient=0.45,
                                                 residual_stomatal_conductance=3.0):
    """Calcualtes the bulk surface conductance of a leaf layer.

    Args:
        incident_par (float): [W m-2ground] incident photosynthetically active radiation above the canopy
        maximum_stomatal_conductance (float): [m h-1] maximum stomatal conductance
        stomatal_sensibility (float): [-] stomatal closure fraction due to water stress
        upper_cumulative_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf area index above the considered
            layer
        lower_cumulative_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf area index below the considered
            layer
        absorbed_par_50 (float): [W m-2leaf] an empirical parameter to regulate the shape of stomatal conductance
            response to absorbed PAR
        global_extinction_coefficient (float): [m2ground m-2leaf] extinction coefficient of the entire PAR wavelength
            into the canopy
        residual_stomatal_conductance (float): [m h-1] residual (minimum) stomatal conductance

    Returns:
        (float): [m h-1] bulk surface conductance of the leaf layer

    """
    scaling_factor = 1.0 / global_extinction_coefficient * log(
        (global_extinction_coefficient * incident_par *
         exp(-global_extinction_coefficient * upper_cumulative_leaf_area_index) + absorbed_par_50) /
        (global_extinction_coefficient * incident_par *
         exp(-global_extinction_coefficient * lower_cumulative_leaf_area_index) + absorbed_par_50))

    return residual_stomatal_conductance * (
            1 - scaling_factor) + maximum_stomatal_conductance * stomatal_sensibility * scaling_factor


def calc_leaf_layer_surface_resistance_to_vapor(incident_par, maximum_stomatal_conductance, stomatal_sensibility,
                                                upper_cumulative_leaf_area_index, lower_cumulative_leaf_area_index,
                                                absorbed_par_50=105, global_extinction_coefficient=0.45,
                                                residual_stomatal_conductance=3.0):
    """Calcualtes the bulk surface resistance of a leaf layer.

    Args:
        incident_par (float): [W m-2ground] incident photosynthetically active radiation above the canopy
        maximum_stomatal_conductance (float): [m h-1] maximum stomatal conductance
        stomatal_sensibility (float): [-] stomatal closure fraction due to water stress
        upper_cumulative_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf area index above the considered
            layer
        lower_cumulative_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf area index below the considered
            layer
        absorbed_par_50 (float): [W m-2leaf] an empirical parameter to regulate the shape of stomatal conductance
            response to absorbed PAR
        global_extinction_coefficient (float): [m2ground m-2leaf] extinction coefficient of PAR into the canopy
        residual_stomatal_conductance (float): [m h-1] residual (minimum) stomatal conductance

    Returns:
        (float): [h m-1] bulk surface resistance of the leaf layer

    """
    surface_conductance = calc_leaf_layer_surface_conductance_to_vapor(incident_par,
                                                                       maximum_stomatal_conductance,
                                                                       stomatal_sensibility,
                                                                       upper_cumulative_leaf_area_index,
                                                                       lower_cumulative_leaf_area_index,
                                                                       absorbed_par_50,
                                                                       global_extinction_coefficient,
                                                                       residual_stomatal_conductance)

    return 1.0 / max(1.e-6, surface_conductance)


def calc_component_net_radiation(net_shortwave_radiation, net_longwave_radiation, is_soil):
    """Calculates the net radiation flux density of a canopy component per unit ground surface area.

    Args:
        net_shortwave_radiation (float): [W m-2ground] net shortwave flux density radiation of the component
        net_longwave_radiation (float): [W m-2ground] net longwave flux density radiation of the component
        is_soil (bool): `True` if the calculation is performed for the soil component, otherwise `False`

    Returns:
        (float): [W m-2ground] net radiation flux density of the canopy component per unit ground surface area

    """
    if not is_soil:
        return net_shortwave_radiation + net_longwave_radiation
    else:
        net_above_ground_radiation = net_shortwave_radiation + net_longwave_radiation
        return net_above_ground_radiation - calc_soil_heat_flux(net_above_ground_radiation, net_shortwave_radiation > 0)


def calc_soil_heat_flux(net_above_ground_radiation, is_diurnal):
    """Calculates the net heat flux density into the soil substrates.

    Args:
        net_above_ground_radiation (float): [W m-2ground] the net above-ground radiation
        is_diurnal (bool): `True` during the daylight hours, otherwise `False`

    Returns:
        (float): [W m-2ground] the net heat flux density into the soil substrates

    """
    return 0.1 * net_above_ground_radiation if is_diurnal else 0.5 * net_above_ground_radiation


def calc_canopy_aerodynamic_resistance(wind_speed, canopy_height, reference_height, von_karman_cst):
    """Calculates air resistance to water vapor transfer between the source height and the reference height.

    Args:
        wind_speed (float): [m h-1] wind speed at measurement height
        canopy_height (float): [m] average height of the canopy
        reference_height (float): [m] height at which wind speed in measured
        von_karman_cst (float): [-] von Karman constant

    Returns:
        (float): [h m-1] air resistance to water vapor transfer between the source height and the reference height

    """
    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)
    zero_displacement_height_ = calc_zero_displacement_height(canopy_height)
    canopy_roughness_length_for_momentum_ = calc_canopy_roughness_length_for_momentum(canopy_height)
    canopy_roughness_length_for_heat_transfer_ = calc_canopy_roughness_length_for_heat_transfer(canopy_height)

    return 1.0 / (wind_speed * von_karman_cst ** 2) * (
            log((reference_height - zero_displacement_height_) / canopy_roughness_length_for_momentum_) *
            log((reference_height - zero_displacement_height_) / canopy_roughness_length_for_heat_transfer_))


def calc_vapor_pressure_slope(air_temperature):
    """Calculates the slope of vapor pressure curve at a given air temperature.

    Args:
        air_temperature (float): [°C] air temperature

    Returns:
        (float): [kPa °K-1] the slope of vapor pressure curve at the given air temperature

    """
    return 4098 * (0.6108 * exp((17.27 * air_temperature) / (air_temperature + 237.3))) / (
            (air_temperature + 237.3) ** 2)


def calc_penman_evaporative_energy(canopy_aerodynamic_resistance, canopy_net_radiation, vapor_pressure_slope,
                                   vapor_pressure_deficit, psychrometric_cst, air_density, air_specific_heat_capacity):
    """Calculates the evapotranspiration energy flux density according to Penman's formula

    Args:
        canopy_aerodynamic_resistance (float): [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        canopy_net_radiation (float): [W m-2ground] net radiation flux density of the entire canopy per unit ground
            surface area
        vapor_pressure_slope (float): [kPa K-1] the slope of vapor pressure curve at a given air temperature
        vapor_pressure_deficit (float): [kPa] vapor pressure deficit at the reference height
        psychrometric_cst (float): [kPa K-1] the psychrometric constant
        air_density (float): [g m-3] density of dry air
        air_specific_heat_capacity (float): [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        (float): [W m-2ground] the evapotranspiration energy flux density according to Penman's formula

    """
    # canopy_net_radiation = max(0.0, canopy_net_radiation)  # to avoid calculating negative transpiration fluxes
    return (vapor_pressure_slope * canopy_net_radiation + (
            air_density * air_specific_heat_capacity * vapor_pressure_deficit) / canopy_aerodynamic_resistance) / (
                   vapor_pressure_slope + psychrometric_cst)


def calc_canopy_lumped_aerodynamic_resistance(canopy_aerodynamic_resistance, vapor_pressure_slope, psychrometric_cst):
    """Calculates Lhomme's lumped aerodynamic resistance (R0).

    Args:
        canopy_aerodynamic_resistance (float): [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        vapor_pressure_slope (float): [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_cst (float): [kPa K-1] the psychrometric constant

    Returns:
        (float): [h m-1] Lhomme's lumped aerodynamic resistance (R0)

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 7

    """
    return (1 + vapor_pressure_slope / psychrometric_cst) * canopy_aerodynamic_resistance


def calc_component_composed_resistance(surface_resistance, boundary_layer_resistance, vapor_pressure_slope,
                                       psychrometric_cst, stomatal_density_factor):
    """Calculates Lhomme's lumped boundary and surface resistance (Ri) for a given canopy component.

    Args:
        surface_resistance (float): [h m-1] surface resistance of the component
        boundary_layer_resistance (float): [h m-1] boundary layer resistance of the component
        vapor_pressure_slope (float): [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_cst (float): [kPa K-1] the psychrometric constant
        stomatal_density_factor (int): [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    Returns:
        (float): [h m-1] Lhomme's lumped boundary and surface resistance (Ri) for the given canopy component

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 8

    """
    return surface_resistance + (
            stomatal_density_factor + vapor_pressure_slope / psychrometric_cst) * boundary_layer_resistance


def calc_component_composed_conductance(composed_boundary_and_surface_resistance,
                                        sum_composed_boundary_and_surface_conductances,
                                        canopy_lumped_aerodynamic_resistance):
    """Calculates Lhomme's composed boundary and surface conductance (Pi) of a given canopy component

    Args:
        composed_boundary_and_surface_resistance (float): [h m-1] Lhomme's lumped boundary and surface resistance for
            the given component
        sum_composed_boundary_and_surface_conductances (float): [m h-1] the sum of Lhomme's lumped boundary and surface
            resistance over all canopy components
        canopy_lumped_aerodynamic_resistance (float): [h m-1] Lhomme's lumped aerodynamic resistance

    Returns:
        (float): [m h-1] Lhomme's composed boundary and surface conductance (Pi) of the given canopy component

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 13


    """
    return 1.0 / (composed_boundary_and_surface_resistance * (
            1.0 + canopy_lumped_aerodynamic_resistance * sum_composed_boundary_and_surface_conductances))


def calc_penman_monteith_evaporative_energy(components_indices, canopy_lumped_aerodynamic_resistance,
                                            penman_evaporative_energy, composed_boundary_and_surface_conductances,
                                            net_radiation_fluxes, boundary_layer_resistances, vapor_pressure_slope,
                                            psychrometric_cst):
    """Calculates canopy evapotrative energy.

    Args:
        components_indices (list): [-] list of dictionary keys of all components of the canopy
        canopy_lumped_aerodynamic_resistance (float): [h m-1] Lhomme's lumped aerodynamic resistance (R0)
        penman_evaporative_energy (float): [W m-2ground] the evapotranspiration energy flux density according to
            Penman's formula
        composed_boundary_and_surface_conductances (dict): [m h-1] Lhomme's composed boundary and surface conductance
            (Pi) of all canopy components
        net_radiation_fluxes (dict): [W m-2ground] net radiation flux densities of all canopy components
        boundary_layer_resistances (dict): [h m-1] boundary layer reistance to heat transfer of all canopy components
        vapor_pressure_slope (float): [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_cst (float): [kPa K-1] the psychrometric constant

    Returns:
        (float): [W m-2ground] canopy evaporative energy

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
        sum_p_a_r += \
            composed_boundary_and_surface_conductances[i] * net_radiation_fluxes[i] * boundary_layer_resistances[i]
    return canopy_lumped_aerodynamic_resistance * penman_evaporative_energy * sum_p + (
            vapor_pressure_slope / psychrometric_cst) * sum_p_a_r


def calc_component_evaporative_energy(net_radiation, boundary_layer_resistance, lumped_boundary_and_surface_resistance,
                                      canopy_lumped_aerodynamic_resistance, penman_evaporative_energy,
                                      penman_monteith_evaporative_energy, vapor_pressure_slope, psychrometric_cst):
    """Calculates the evaporative energy for a given canopy component.

    Args:
        net_radiation (float): [W m-2ground] net radiation flux density of the given component
        boundary_layer_resistance (float): [h m-1] component boundary layer resistance to heat transfer
        lumped_boundary_and_surface_resistance (float): [h m-1] Lhomme's lumped boundary and surface resistance (Ri)
            for the given component
        canopy_lumped_aerodynamic_resistance (float): [h m-1] Lhomme's lumped aerodynamic resistance (R0)
        penman_evaporative_energy (float): [W m-2ground] the evapotranspiration energy flux density according to
            Penman's formula
        penman_monteith_evaporative_energy (float): [W m-2ground] canopy evaporative energy
        vapor_pressure_slope (float): [kPa K-1] the slope of vapor pressure curve at a given air temperature
        psychrometric_cst (float): [kPa K-1] the psychrometric constant

    Returns:
        (float): [W m-2ground] the evaporative energy for the given canopy component

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 11

    """
    energy_driven_evaporation = canopy_lumped_aerodynamic_resistance * (
            penman_evaporative_energy - penman_monteith_evaporative_energy)
    radiation_driven_evaporation = boundary_layer_resistance * net_radiation * (
            vapor_pressure_slope / psychrometric_cst)
    return 1.0 / lumped_boundary_and_surface_resistance * (energy_driven_evaporation + radiation_driven_evaporation)


def calc_canopy_temperature(air_temperature, canopy_aerodynamic_resistance, canopy_net_radiation,
                            penman_monteith_evaporative_energy, air_density, air_specific_heat_capacity):
    """Calculates air temperature at source height.

    Args:
        air_temperature (float): [K] air temperature at measurement height
        canopy_aerodynamic_resistance (float): [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        canopy_net_radiation (float): [W m-2ground] canopy net radiation
        penman_monteith_evaporative_energy (float): [W m-2ground] canopy evaporative energy
        air_density (float): [g m-3] density of dry air
        air_specific_heat_capacity (float): [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        (float): [K] air temperature at source height

    """
    return air_temperature + (canopy_aerodynamic_resistance / (air_density * air_specific_heat_capacity)) * (
            canopy_net_radiation - penman_monteith_evaporative_energy)


def calc_component_temperature(canopy_temperature, boundary_layer_resistance, component_net_radiation,
                               component_evaporative_energy, air_density, air_specific_heat_capacity):
    """Calculates air temperature at source height.

    Args:
        canopy_temperature (float): [K] air temperature at source height
        boundary_layer_resistance (float): [h m-1] resistance to heat transfer at the given component
        component_net_radiation (float): [W m-2ground] canopy net radiation
        component_evaporative_energy (float): [W m-2ground] canopy evaporative energy
        air_density (float): [g m-3] density of dry air
        air_specific_heat_capacity (float): [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        (float): [K] air temperature at source height

    """
    return canopy_temperature + (boundary_layer_resistance / (air_density * air_specific_heat_capacity)) * (
            component_net_radiation - component_evaporative_energy)


def calc_temperature_step(previous_value, actual_value, step_fraction=0.5):
    """Calculates the temperature step value between two consecutive energy balance calculations.

    Args:
        previous_value (float): [K] previous calculated temperature
        actual_value (float): [K] actual calculated temperature
        step_fraction (float): [-] fraction of the entire step (`actual_value - previous_value`) to be used

    Returns:
        (float): [K] the temperature step value between two consecutive energy balance calculations

    """
    return step_fraction * (actual_value - previous_value)


def calc_component_net_longwave_radiation(air_temperature, lower_leaf_area_index, atmoshperic_emissivity,
                                          extinction_coef, stefan_boltzman_cst):
    """Calculates net long wave radiation of a canopy component.

    Args:
        air_temperature (float): [K] air temperature
        lower_leaf_area_index (float): [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        atmoshperic_emissivity (float): [-] atmospheric emissivity to longwave radiation
        extinction_coef (float): [m2ground m-2leaf] extinction coefficient of longwave radiation, assumed equal to
            the extinction coefficient of the diffuse photosynthetically active radiation for black leaves
        stefan_boltzman_cst (float):

    Returns:

    """
    scaling_factor = 1.0 / extinction_coef * exp(-extinction_coef * lower_leaf_area_index)
    return - (1 - atmoshperic_emissivity) * stefan_boltzman_cst * air_temperature ** 4 * scaling_factor
