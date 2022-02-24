from crop_energy_balance.formalisms.config import PRECISION


def calc_boundary_layer_resistance(forced_convection_resistance: float,
                                   free_convection_resistance: float) -> float:
    """Calculates total boundary resistance under both forced and free convection conditions.

    Args:
        forced_convection_resistance: [h m-1] boundary resistance under forced convection
        free_convection_resistance: [h m-1] boundary resistance under free convection

    Returns:
        [h m-1] boundary resistance under both forced and free convection conditions
    """

    return 1. / (1. / max(PRECISION, forced_convection_resistance) + 1. / max(PRECISION, free_convection_resistance))


def calc_composed_resistance(surface_resistance: float,
                             boundary_layer_resistance: float,
                             vapor_pressure_slope: float,
                             psychrometric_constant: float,
                             stomatal_density_factor: int) -> float:
    """Calculates Lhomme's lumped boundary and surface resistance (Ri) for a given canopy component.

    Args:
        surface_resistance: [h m-1] surface resistance to water vapor transfer
        boundary_layer_resistance: [h m-1] boundary layer resistance to heat transfer
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


def calc_composed_conductance(composed_boundary_and_surface_resistance: float,
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


def calc_evaporative_energy(available_energy: float,
                            boundary_layer_resistance: float,
                            lumped_boundary_and_surface_resistance: float,
                            canopy_lumped_aerodynamic_resistance: float,
                            penman_evaporative_energy: float,
                            penman_monteith_evaporative_energy: float,
                            vapor_pressure_slope: float,
                            psychrometric_constant: float) -> float:
    """Calculates the evaporative energy for a given canopy component.

    Args:
        available_energy: [W m-2ground] available energy flux density of the given component
        boundary_layer_resistance: [h m-1] boundary layer resistance to heat transfer
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
    radiation_driven_evaporation = boundary_layer_resistance * available_energy * (
            vapor_pressure_slope / psychrometric_constant)
    return 1.0 / lumped_boundary_and_surface_resistance * (energy_driven_evaporation + radiation_driven_evaporation)


def calc_temperature(canopy_temperature: float,
                     boundary_layer_resistance: float,
                     available_energy: float,
                     evaporative_energy: float,
                     air_density: float,
                     air_specific_heat_capacity: float) -> float:
    """Calculates air temperature at source height.

    Args:
        canopy_temperature: [K] air temperature at source height
        boundary_layer_resistance: [h m-1] boundary layer resistance to heat transfer
        available_energy: [W m-2ground] available energy to the component
        evaporative_energy: [W m-2ground] canopy evaporative energy
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [K] air temperature at source height
    """
    return canopy_temperature + (boundary_layer_resistance / (air_density * air_specific_heat_capacity)) * (
            available_energy - evaporative_energy)


def calc_available_energy(net_shortwave_radiation: float,
                          net_longwave_radiation: float,
                          soil_heat_flux: float) -> float:
    """Calculates the available energy flux density of a canopy component per unit ground surface area.

    Args:
        net_shortwave_radiation: [W m-2ground] net shortwave flux density radiation of the component
        net_longwave_radiation: [W m-2ground] net longwave flux density radiation of the component
        soil_heat_flux: [W m-2ground] the net heat flux density into the soil

    Returns:
        [W m-2ground] available energy flux density of the component per unit ground surface area
    """
    return net_shortwave_radiation + net_longwave_radiation - soil_heat_flux
