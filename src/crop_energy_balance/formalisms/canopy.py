from math import log


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
    wind_speed = max(2400., wind_speed)
    canopy_height = max(0.1, canopy_height)

    return (von_karman_constant ** 2 * wind_speed * (canopy_height - zero_displacement_height)) / (
        log((measurement_height - zero_displacement_height) / canopy_roughness_length_for_momentum))


def calc_net_longwave_radiation(air_temperature: float,
                                atmospheric_emissivity: float,
                                stefan_boltzman_constant: float) -> float:
    """Calculates net long wave radiation at the top of the canopy.

    Args:
        air_temperature: [K] air temperature
        atmospheric_emissivity: [-] atmospheric emissivity to longwave radiation
        stefan_boltzman_constant: [W m-2 K-4] Stefan-Boltzmann constant

    Returns:
        [W m-2ground] net long wave radiation at the top of the canopy

    References:
        Leuning et al. 1995
            Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies.
            Plant, Cell and Environment 18, 1183 - 1200.
    """
    return (1 - atmospheric_emissivity) * stefan_boltzman_constant * air_temperature ** 4


def calc_canopy_aerodynamic_resistance(wind_speed: float,
                                       canopy_height: float,
                                       reference_height: float,
                                       von_karman_constant: float) -> float:
    """Calculates air resistance to water vapor transfer between the source height and the reference height.

    Args:
        wind_speed: [m h-1] wind speed at measurement height
        canopy_height: [m] average height of the canopy
        reference_height: [m] height at which wind speed in measured
        von_karman_constant: [-] von Karman constant

    Returns:
        [h m-1] air resistance to water vapor transfer between the source height and the reference height
    """
    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)
    zero_displacement_height_ = calc_zero_displacement_height(canopy_height)
    canopy_roughness_length_for_momentum_ = calc_canopy_roughness_length_for_momentum(canopy_height)
    canopy_roughness_length_for_heat_transfer_ = calc_canopy_roughness_length_for_heat_transfer(canopy_height)

    return 1.0 / (wind_speed * von_karman_constant ** 2) * (
            log((reference_height - zero_displacement_height_) / canopy_roughness_length_for_momentum_) *
            log((reference_height - zero_displacement_height_) / canopy_roughness_length_for_heat_transfer_))


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


def calc_penman_monteith_evaporative_energy(canopy_lumped_aerodynamic_resistance: float,
                                            penman_evaporative_energy: float,
                                            composed_boundary_and_surface_conductances: list,
                                            net_radiation_fluxes: list,
                                            boundary_layer_resistances: list,
                                            vapor_pressure_slope: float,
                                            psychrometric_constant: float) -> float:
    """Calculates canopy evaporative energy.

    Args:
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

    Notes:
        All lists of :args:`composed_boundary_and_surface_conductances`, :args:`net_radiation_fluxes` and
            :args:`boundary_layer_resistances` must be equally ordered, that is, their values must refer to the same crop
            components in the same order. Violating this condition will lead serious to simulation errors.
    """
    sum_p = 0.0
    sum_p_a_r = 0.0
    for p, a, r in zip(composed_boundary_and_surface_conductances, net_radiation_fluxes, boundary_layer_resistances):
        sum_p += p
        sum_p_a_r += (p * a * r)
    return canopy_lumped_aerodynamic_resistance * penman_evaporative_energy * sum_p + (
            vapor_pressure_slope / psychrometric_constant) * sum_p_a_r


def calc_temperature(air_temperature: float,
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
