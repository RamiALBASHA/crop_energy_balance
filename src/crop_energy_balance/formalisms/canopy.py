from math import log, atan, pi

from crop_energy_balance.formalisms.config import PRECISION
from crop_energy_balance.formalisms.weather import convert_celsius_to_kelvin


def calc_zero_displacement_height(canopy_height: float, leaf_area_index: float, drag_coefficient: float) -> float:
    """Calculates zero displacement height of the canopy.

    Args:
        canopy_height: [m] average height of the canopy
        leaf_area_index: [m2leaf m-2ground] total leaf area index
        drag_coefficient: [m2ground m-2leaf] drag coefficient

    Returns:
        [m] zero displacement height of water vapor between the crop and the atmosphere

    References:
        Choudhury, B., Monteith, J., 1988.
            A four-layer model for the heat budget ofhomogeneous land surfaces.
            Q. J. Roy. Meteor. Soc. 114, 373–398.

    """
    return 1.1 * canopy_height * log(1.0 + (drag_coefficient * leaf_area_index) ** 0.25)


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


def calc_roughness_length_for_heat_transfer(roughness_length_for_momentum_transfer: float,
                                            ratio_heat_to_momentum_roughness_lengths: float) -> float:
    """Calculates roughness length for heat and water vapor transfer of the canopy.

    Args:
        roughness_length_for_momentum_transfer: [m] roughness length for momentum transfer
        ratio_heat_to_momentum_roughness_lengths: [-] Ratio of canopy's heat to momentum roughness lengths

    Returns:
        [m] roughness length for heat transfer

    References:
        Shuttleworth 2007.
            Putting the 'vap' into evaporation.
            Hydrol. Earth Syst. Sci. 11, 210 - 244
            Eq. 11
        Kimball et al. 2015.
            Predicting canopy temperatures and infrared heater energy requirements for warming field plots.
            Climatology and Water Management 107, 129 - 141
            Eq. 7

    Notes:
        Indicative values for 'ratio_heat_to_momentum_roughness_lengths' are:
            * 1/10 for reference grass crop (Shuttleworth, 2007)
            * 1/7.4 for wheat (Kimball et al., 2015)

    """
    return ratio_heat_to_momentum_roughness_lengths * roughness_length_for_momentum_transfer


def calc_wind_speed_at_canopy_height(wind_speed: float,
                                     canopy_height: float,
                                     measurement_height: float,
                                     zero_displacement_height: float,
                                     roughness_length_for_momentum: float) -> float:
    """Calculates hourly wind speed values at canopy height.

    Args:
        wind_speed: [m h-1] wind speed at measurement height
        canopy_height: [m] height of the canopy
        measurement_height: [m] height at which meteorological measurements are made
        zero_displacement_height: [m] zero displacement height
        roughness_length_for_momentum: [m] roughness length for momentum transfer

    Returns:
        [m h-1] wind speed at canopy height
    """

    canopy_height = max(0.1, canopy_height)

    return max(PRECISION,
               wind_speed * log((canopy_height - zero_displacement_height) / roughness_length_for_momentum) / log(
                   (measurement_height - zero_displacement_height) / roughness_length_for_momentum))


def calc_net_longwave_radiation(air_temperature: float,
                                air_vapor_pressure: float,
                                atmospheric_emissivity: float,
                                stefan_boltzmann_constant: float) -> float:
    """Calculates net long wave radiation at the top of the canopy.

    Args:
        air_temperature: [K] air temperature
        air_vapor_pressure: [kPa] actual vapor pressure
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
    gross_radiation_loss = - stefan_boltzmann_constant * air_temperature ** 4
    correction_humidity = 0.34 - 0.14 * air_vapor_pressure ** 0.5
    correction_sky_cover = (1. - atmospheric_emissivity)
    return gross_radiation_loss * correction_humidity * correction_sky_cover


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
                                richardon_threshold_free_convection: float,
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
        richardon_threshold_free_convection: [-] Richardson number threshold below which free convection is assumed.
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
    if richardson_number < richardon_threshold_free_convection:
        # ----------------------
        # strongly unstable, free convection dominates
        # (Webber et al. 2016, eq. 11 is replaced by Kimball et al. 2015 eqs. 28 & 29)
        # Kimball et al. (2015) used the equation for free (also called natural) convection from ASHAE (page 3.12,
        # Table 5, eq. 11). They replaced the "n" parameter value from 1.52 used for "large plates" to 5 for wheat.
        # Note that the value of h (Heat Transfer Coefficient) given in ASHAE is in W m-2 K-1.
        # ----------------------
        n = 5.
        temperature_difference = max(convert_celsius_to_kelvin(0.1), abs(canopy_temperature - air_temperature))
        ra = air_density * air_specific_heat_capacity / (n * temperature_difference ** (1 / 3.))
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


def calc_friction_velocity(wind_speed: float,
                           measurement_height: float,
                           zero_displacement_height: float,
                           roughness_length_for_momentum: float,
                           stability_correction_for_momentum: float,
                           von_karman_constant: float):
    """Calculates the friction velocity in the atmospheric boundary layer.

    Args:
        wind_speed: [m h-1] wind speed at measurement height
        measurement_height: [m] height at which meteorological measurements are made
        zero_displacement_height: [m] zero plane displacement height
        roughness_length_for_momentum: [m] roughness length for momentum transfer
        stability_correction_for_momentum: [-] stability correction factor for momentum transfer
        von_karman_constant: [-] von Karman constant

    Returns:
        [m h-1]: friction velocity
    """
    return von_karman_constant * wind_speed / (
            log((measurement_height - zero_displacement_height) / roughness_length_for_momentum)
            - stability_correction_for_momentum)


def calc_monin_obukhov_length(surface_temperature: float,
                              sensible_heat_flux: float,
                              friction_velocity: float,
                              air_density: float,
                              air_specific_heat_capacity: float,
                              gravitational_acceleration: float,
                              von_karman_constant: float) -> float:
    """Calculates the Monin-Obukhov length.

    Args:
        surface_temperature: [K] aerodynamic surface temperature
        sensible_heat_flux: [W m-2ground] sensible heat flux
        friction_velocity: [m h-1] friction velocity
        air_density: [g m-3] air density
        air_specific_heat_capacity: [W h g-1 K-1] air specific heat capacity
        gravitational_acceleration: [m h-2] gravity acceleration
        von_karman_constant: [-] von Karman constant

    Returns:
        [m] Monin-Obukhov length

    References:
        Monin, A. S. and Obukhov, A. M. 1954.
            Dimensionless characteristics of turbulence in the surface layer.
            Akad. Nauk. SSSR Geofiz. Inst. Tr. 24, 163 - 187.
    """
    shear = friction_velocity ** 3
    buoyancy = von_karman_constant * (gravitational_acceleration / surface_temperature) * (
            sensible_heat_flux / (air_density * air_specific_heat_capacity))
    return - shear / buoyancy


def calc_richardson_number(is_stable: bool,
                           measurement_height: float,
                           zero_displacement_height: float,
                           monin_obukhov_length: float) -> float:
    """Calculates the Richardson's number.

    Args:
        is_stable: first approximation of stability (see Note below)
        measurement_height: [m] height at which meteorological measurements are made
        zero_displacement_height: [m] zero plane displacement height
        monin_obukhov_length: [m] Monin-Obukhov length

    Returns:
        [-] Richardson number

    References:
        Monteith and Unsworth (2013).
            Principles of Environmental Physics (Fourth Edition)
            Academic Press, pp 289 - 320

        Webb (1970).
            Profile relationships: the log-linear range, and extension to strong stability.
            Quarterly Journal of the Royal Meteorological Society 96, 67 - 90.

        Webber et al. (2016)
            Simulating canopy temperature for modelling heat stress in cereals.
            Environmental Modelling and Software 77, 143 - 155


    Note:
        According to Webber et al. (2016) turbulence condition is considered stable if canopy's temperature is lower
            the air's, otherwise unstable
    """
    if is_stable:
        # Webb (1970) in Monteith and Unsworth (2013)
        richardson = (measurement_height - zero_displacement_height) / (
                monin_obukhov_length + 5 * (measurement_height - zero_displacement_height))
    else:
        # Monteith and Unsworth (2013)
        richardson = (measurement_height - zero_displacement_height) / monin_obukhov_length

    return richardson


def calc_stability_correction_functions(monin_obukhov_length: float,
                                        richardson_number: float,
                                        measurement_height: float,
                                        zero_displacement_height: float,
                                        richardon_threshold_free_convection: float) -> (float,):
    """Calculates the stability correction functions for momentum and heat transfer.

    Args:
        monin_obukhov_length: [m] Monin-Obukhov length
        richardson_number: [-] Richardson length
        measurement_height: [m] measurement height for wind and temperature measurement (assumed equal)
        zero_displacement_height: [m] zero plane displacement height
        richardon_threshold_free_convection: [-] Richardson number threshold below which free convection is assumed

    Returns:
        [-] correction function for momentum transfer
        [-] correction function for heat transfer
    """

    if richardson_number < richardon_threshold_free_convection:
        # ----------------------
        # strongly unstable, free convection dominates
        # ----------------------
        correction_for_heat = 0
        correction_for_momentum = correction_for_heat
    elif richardson_number < -0.01:
        # ----------------------
        # unstable
        # ----------------------
        # Colaizzi et al. 2004, eq. 12)
        x = (1.0 - 16.0 * (measurement_height - zero_displacement_height) / monin_obukhov_length) ** 0.25
        correction_for_heat = 2.0 * log((1 + x ** 2) / 2)  # (Liu et al. 2007, eq. 13)
        correction_for_momentum = 2.0 * log((1 + x) / 2) + log((1 + x ** 2) / 2) - 2 * atan(x) + pi / 2.
    elif richardson_number < 0.2:
        # ----------------------
        # stable - technically (Thom, 1975)
        # ----------------------
        correction_for_heat = -5.0 * (
                measurement_height - zero_displacement_height) / monin_obukhov_length  # (Liu et al. 2007, eq. 11)
        correction_for_momentum = correction_for_heat  # (Liu et al. 2007, eq. 10)
    else:
        # ----------------------
        # strongly stable
        # Note: The concepts underlying similarity theory become invalid according to Mahrt (2010)
        # in Monteith and Unsworth (2013)
        # ----------------------
        correction_for_heat = 0
        correction_for_momentum = correction_for_heat

    return correction_for_momentum, correction_for_heat
