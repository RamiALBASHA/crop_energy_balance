from math import log

from crop_energy_balance.formalisms.config import PRECISION


def calc_zero_displacement_height(canopy_height: float) -> float:
    """Calculates zero displacement height of the canopy.

    Args:
        canopy_height: [m] average height of the canopy

    Returns:
        [m] zero displacement height of water vapor between the crop and the atmosphere
    """
    return 0.67 * canopy_height


def calc_roughness_length_for_momentum_transfer(canopy_height: float) -> float:
    """Calculates roughness length for momentum transfer of the canopy.

    Args:
        canopy_height: [m] average height of the canopy

    Returns:
        [m] roughness length for momentum

    References:
        Shuttleworth 2007.
            Putting the 'vap' into evaporation.
            Hydrol. Earth Syst. Sci. 11, 210 - 244
            Eq. 11

    Notes:
        This formula holds for the reference grass crop.
    """
    return 0.123 * canopy_height


def calc_roughness_length_for_heat_transfer(canopy_height: float) -> float:
    """Calculates roughness length for heat and water vapor transfer of the canopy.

    Args:
        canopy_height: [m] average height of the canopy

    Returns:
        [m] roughness length for heat transfer

    References:
        Shuttleworth 2007.
            Putting the 'vap' into evaporation.
            Hydrol. Earth Syst. Sci. 11, 210 - 244
            Eq. 11

    Notes:
        This formula holds for the reference grass crop.
    """
    return calc_roughness_length_for_momentum_transfer(canopy_height=canopy_height) / 10.


def calc_wind_speed_at_canopy_height(wind_speed: float,
                                     canopy_height: float,
                                     measurement_height: float) -> float:
    """Calculates hourly wind speed values at canopy height.

    Args:
        wind_speed: [m h-1] wind speed at measurement height
        canopy_height: [m] height of the canopy
        measurement_height: [m] height at which meteorological measurements are made

    Returns:
        [m h-1] wind speed at canopy height
    """

    wind_speed = max(2400.0, wind_speed)
    canopy_height = max(0.1, canopy_height)
    d = calc_zero_displacement_height(canopy_height)
    z0u = calc_roughness_length_for_momentum_transfer(canopy_height)

    return max(PRECISION, wind_speed * log((canopy_height - d) / z0u) / log((measurement_height - d) / z0u))


def calc_net_longwave_radiation(air_temperature: float,
                                air_vapor_pressure: float,
                                canopy_temperature: float,
                                atmospheric_emissivity: float,
                                stefan_boltzmann_constant: float) -> float:
    """Calculates net long wave radiation at the top of the canopy.

    Args:
        air_temperature: [K] air temperature
        air_vapor_pressure: [kPa] actual vapor pressure
        canopy_temperature: [K] canopy (source) temperature
        atmospheric_emissivity: [-] atmospheric emissivity to longwave radiation
        stefan_boltzmann_constant: [W m-2 K-4] Stefan-Boltzmann constant

    Returns:
        [W m-2ground] net long wave radiation at the top of the canopy

    References:
        Allen et al. 1998
            Crop Evapotranspiration – Guide-lines for Computing Crop Water Requirements
            Paper 56. Food and Agricultural Organization of the United Nations.
        Leuning et al. 1995
            Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies.
            Plant, Cell and Environment 18, 1183 - 1200.
    """
    gross_radiation_loss = (atmospheric_emissivity * stefan_boltzmann_constant * air_temperature ** 4 -
                            stefan_boltzmann_constant * canopy_temperature ** 4)
    reduction_factor = 0.34 - 0.14 * air_vapor_pressure ** 0.5
    return gross_radiation_loss * reduction_factor


def calc_sensible_heat_flux(source_temperature: float,
                            air_temperature: float,
                            aerodynamic_resistance: float,
                            air_density: float,
                            air_specific_heat_capacity: float) -> float:
    """Calculates the sensible heat flux of the canopy.

    Args:
        source_temperature: [K] air temperature at source height
        air_temperature: [K] air temperature at reference height
        aerodynamic_resistance: [h m-1] air resistance to heat transfer between the source height and the reference
            height
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [W m-2ground]: Sensible heat flux of the canopy.
    """
    temperature_difference = source_temperature - air_temperature
    return air_density * air_specific_heat_capacity * temperature_difference / aerodynamic_resistance


def calc_aerodynamic_resistance(richardson_number: float,
                                friction_velocity: float,
                                measurement_height: float,
                                zero_displacement_height: float,
                                roughness_length_for_heat: float,
                                stability_correction_for_heat: float,
                                canopy_temperature: float,
                                air_temperature: float,
                                von_karman_constant: float,
                                air_density: float,
                                air_specific_heat_capacity: float):
    """Calculates the resistance to the transfer of momentum, water vapor and heat between the source height and the
    reference height.

    Args:
        friction_velocity: [m h-1] friction velocity
        richardson_number: [-] Richardson number (indicates stability condition)
        measurement_height: [m] height at which meteorological measurements are made
        zero_displacement_height: [m] zero displacement height
        roughness_length_for_heat: [m] roughness length for heat and water vapor transfer
        stability_correction_for_heat: [-] stability correction factor for heat transfer
        von_karman_constant: [-] von Karman constant
        air_temperature: [K] air temperature at reference height
        canopy_temperature: [K] canopy (source) temperature
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [h m-1] air resistance to water vapor transfer between the source height and the reference height

    Notes:
        For richardson_number < -0.8, the aerodynamic resistance is calculated for free convection conditions following
            Kimball et al. (2015) reported by Webber et al. (2016). The characteristic surface factor
            (n in Webber et al. 2016) was set to 5.0 as suggested by Kimball et al. (2015, Eq. 28) for open wheat field.
            Yet, originally, n = 1.52 for large plates (ASHRAE, 2001).

    References:
        American Society of Heating, Refrigerating, and Air-Conditioning Engineers (2001).
            ASHRAE Fundamentals Handbook (SI). ASHRAE, New York.
        Kimball et al. (2015)
            Predicting Canopy Temperatures and Infrared Heater Energy Requirements for Warming Field Plots.
            Agronomy Journal 107, 129 - 141.
        Webber et al. (2016)
            Simulating canopy temperature for modelling heat stress in cereals.
            Environmental Modelling and Software 77, 143 - 155

    """
    if richardson_number < -0.8:
        # ----------------------
        # strongly unstable, free convection dominates (Kimball et al. 2015)
        # (Webber et al. 2016, eq. 11)
        # ----------------------
        n = 5.  # surface characteristic shape factor (ASHRAE, 2001, Eq. 11 in Table 5)
        ra = air_density * air_specific_heat_capacity / (n * abs((canopy_temperature - air_temperature)) ** 0.25)
    else:
        # ----------------------
        # unstable and stable conditions (the general case)
        # ----------------------
        ra = 1.0 / (von_karman_constant * friction_velocity) * (
            (log((measurement_height - zero_displacement_height) / roughness_length_for_heat) -
             stability_correction_for_heat))
    return ra


def calc_penman_evaporative_energy(canopy_aerodynamic_resistance: float,
                                   canopy_available_energy: float,
                                   vapor_pressure_slope: float,
                                   vapor_pressure_deficit: float,
                                   psychrometric_constant: float,
                                   air_density: float,
                                   air_specific_heat_capacity: float) -> float:
    """Calculates the evapotranspiration energy flux density according to Penman's formula.

    Args:
        canopy_aerodynamic_resistance: [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        canopy_available_energy: [W m-2ground] available energy flux density of the entire canopy
        vapor_pressure_slope: [kPa K-1] the slope of vapor pressure curve at a given air temperature
        vapor_pressure_deficit: [kPa] vapor pressure deficit at the reference height
        psychrometric_constant: [kPa K-1] the psychrometric constant
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [W m-2ground] the evapotranspiration energy flux density according to Penman's formula

    References:
        Lhomme et al. 2013
            Evaporation from multi-component canopies: Generalized formulations.
            Journal of Hydrology 486, 315 - 320.
            Eq. 10
    """
    return (vapor_pressure_slope * canopy_available_energy + (
            air_density * air_specific_heat_capacity * vapor_pressure_deficit) / canopy_aerodynamic_resistance) / (
                   vapor_pressure_slope + psychrometric_constant)


def calc_lumped_aerodynamic_resistance(canopy_aerodynamic_resistance: float,
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
                                            available_energy_fluxes: list,
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
        available_energy_fluxes: [W m-2ground] available energy flux densities of all canopy components
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
        All lists of :args:`composed_boundary_and_surface_conductances`, :args:`available_energy_fluxes` and
            :args:`boundary_layer_resistances` must be identically ordered, that is, their values must refer to the
            same crop components in the same order. Violating this condition will lead serious to simulation errors.
    """
    sum_p = 0.0
    sum_p_a_r = 0.0
    for p, a, r in zip(composed_boundary_and_surface_conductances, available_energy_fluxes, boundary_layer_resistances):
        sum_p += p
        sum_p_a_r += (p * a * r)
    return canopy_lumped_aerodynamic_resistance * penman_evaporative_energy * sum_p + (
            vapor_pressure_slope / psychrometric_constant) * sum_p_a_r


def calc_temperature(air_temperature: float,
                     canopy_aerodynamic_resistance: float,
                     canopy_available_energy: float,
                     penman_monteith_evaporative_energy: float,
                     air_density: float,
                     air_specific_heat_capacity: float) -> float:
    """Calculates air temperature at source height.

    Args:
        air_temperature: [K] air temperature at measurement height
        canopy_aerodynamic_resistance: [h m-1] air resistance to water vapor transfer between the source height
            and the reference height
        canopy_available_energy: [W m-2ground] canopy available energy
        penman_monteith_evaporative_energy: [W m-2ground] canopy evaporative energy
        air_density: [g m-3] density of dry air
        air_specific_heat_capacity: [W h g-1 K-1] specific heat capacity of the air under a constant pressure

    Returns:
        [K] air temperature at source height
    """
    return air_temperature + (canopy_aerodynamic_resistance / (air_density * air_specific_heat_capacity)) * (
            canopy_available_energy - penman_monteith_evaporative_energy)
