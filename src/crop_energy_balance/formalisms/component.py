from math import exp
from crop_energy_balance.formalisms.soil import calc_soil_heat_flux


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


def calc_component_net_longwave_radiation(air_temperature: float,
                                          lower_cumulative_leaf_area_index: float,
                                          atmospheric_emissivity: float,
                                          extinction_coefficient: float,
                                          stefan_boltzman_constant: float) -> float:
    """Calculates net long wave radiation of a canopy component.

    Args:
        air_temperature: [K] air temperature
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        atmospheric_emissivity: [-] atmospheric emissivity to longwave radiation
        extinction_coefficient: [m2ground m-2leaf] extinction coefficient of longwave radiation, assumed equal to
            the extinction coefficient of the diffuse photosynthetically active radiation for black leaves
        stefan_boltzman_constant: [W m-2 K-4] Stefan-Boltzmann constant

    Returns:
        net long wave radiation of a canopy component
    """
    scaling_factor = 1.0 / extinction_coefficient * exp(-extinction_coefficient * lower_cumulative_leaf_area_index)
    return - (1 - atmospheric_emissivity) * stefan_boltzman_constant * air_temperature ** 4 * scaling_factor


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
