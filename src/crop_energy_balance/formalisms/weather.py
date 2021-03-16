from math import exp, log, atan


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


def calc_atmospheric_emissivity(model: str,
                                air_vapor_pressure: float,
                                air_temperature: float) -> float:
    """Calculates the effective atmospheric emissivity of a clear sky

    Args:
        model: model to be used, one of 'brutsaert_1975' and 'monteith_2013'
        air_vapor_pressure: [kPa] air vapor pressure
        air_temperature: [K] air temperature

    Returns:
        [-] effective atmospheric emissivity of a clear sky

    References:
        Brutsaert, 1975.
            On a Derivable Formula for Long-Wave Radiation From Clear Skies.
            Water Resources Research 11, 742 - 744.
        Monteith and Unsworth, 2013
            Principals of Environmental Physics. Fourth Edition
            pp 72
    """

    if model == 'brutsaert_1975':
        res = 1.24 * (0.1 * air_vapor_pressure / air_temperature) ** (1. / 7.)
    elif model == 'monteith_2013':
        a = 0.10  # kg−1 m−2
        b = 1.2
        c = 0.30  # kg−1 m2
        w = 4.65 * 1.e3 * air_vapor_pressure / air_temperature
        res = 1 - (1 + a * w) * exp(-(b + c * w) ** 0.5)
    else:
        raise ValueError(f'Unknown model name: {model}.')
    return res


def calc_vapor_pressure_slope(temperature: float) -> float:
    """Calculates the slope of air vapor pressure curve at a given air temperature.

    Args:
        temperature: [°C] air temperature

    Returns:
        [kPa K-1] the slope of vapor pressure curve at the given air temperature

    """
    return 4098 * (0.6108 * exp((17.27 * temperature) / (temperature + 237.3))) / ((temperature + 237.3) ** 2)


def convert_kelvin_to_celsius(temperature: float,
                              absolute_zero: float = -273.15) -> float:
    """Converts Kelvin unit to Celsius unit.

    Args:
        temperature: [K] temperature
        absolute_zero: [°C] absolute zero temperature

    Returns:
        [°C] temperature

    """
    return temperature + absolute_zero


def convert_celsius_to_kelvin(temperature: float,
                              absolute_zero: float = -273.15) -> float:
    """Converts Kelvin unit to Celsius unit.

    Args:
        temperature: [°C] temperature
        absolute_zero: [°C] absolute zero temperature

    Returns:
        [K] temperature

    """
    return temperature - absolute_zero


def convert_photosynthetically_active_radiation_into_global_radiation(value: float) -> float:
    """Converts photosynthetically active radiation into global radiation.

    Args:
        value: [W_{PAR} m-2] photosynthetically active radiation

    Returns:
        [W_{global} m-2] global radiation
    """
    return value / 0.48


def convert_global_irradiance_into_photosynthetically_active_radiation(value: float) -> float:
    """Converts global radiation into photosynthetically active radiation.

    Args:
        value: [W_{global} m-2] global radiation

    Returns:
        [W_{PAR} m-2] photosynthetically active radiation
    """
    return value * 0.48


def calc_saturated_air_vapor_pressure(temperature: float) -> float:
    """Compute saturated air vapor pressure.
    Args:
        temperature: [°C] air temperature
    Returns:
        [kPa] saturated air vapor pressure
    """
    return 0.611 * exp(17.27 * temperature / (237.3 + temperature))


def calc_vapor_pressure_deficit(temperature_air: float, temperature_leaf: float, relative_humidity: float) -> float:
    """Computes leaf-to-air vapour pressure deficit.
    Args:
        temperature_air: [°C] air temperature
        temperature_leaf: [°C] leaf temperature
        relative_humidity: [-] air relative humidity (%, between 0 and 100)
    Returns:
        [kPa] leaf-to-air vapour pressure deficit
    """
    es_l = calc_saturated_air_vapor_pressure(temperature_leaf)
    es_a = calc_saturated_air_vapor_pressure(temperature_air)
    ea = es_a * relative_humidity / 100

    return es_l - ea


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
        [-] Richardson length

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


def calc_stability_correction_functions(friction_velocity: float,
                                        sensible_heat: float,
                                        canopy_temperature: float,
                                        measurement_height: float,
                                        zero_displacement_height: float,
                                        air_density: float,
                                        air_specific_heat_capacity: float,
                                        von_karman_constant,
                                        gravitational_acceleration: float) -> (float,):
    """Calculates the stability correction functions for momentum and heat transfer.

    Args:
        friction_velocity: [m h-1] friction velocity
        sensible_heat: [W m-2ground] sensible heat flux
        canopy_temperature: [K] canopy temperature
        measurement_height: [m] measurement height for wind and temperature measurement (assumed equal)
        zero_displacement_height: [m] zero plane displacement height
        air_density: [g m-3] air density
        air_specific_heat_capacity: [W h g-1 K-1] air specific heat capacity
        von_karman_constant: [-] von Karman's constant
        gravitational_acceleration: [m h-2]

    Returns:
        [-] correction function for momentum transfer
        [-] correction function for heat transfer
        [-] Richardson's number
    """

    monin_obukhov_length = calc_monin_obukhov_length(
        surface_temperature=canopy_temperature,
        sensible_heat_flux=sensible_heat,
        friction_velocity=friction_velocity,
        air_density=air_density,
        air_specific_heat_capacity=air_specific_heat_capacity,
        gravitational_acceleration=gravitational_acceleration,
        von_karman_constant=von_karman_constant)

    richardson_number = calc_richardson_number(
        is_stable=sensible_heat < 0,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        monin_obukhov_length=monin_obukhov_length)

    if richardson_number < -0.8:
        # ----------------------
        # strongly unstable, free convection dominates (Kimball et al. 2015)
        # ----------------------
        correction_for_heat = 0
        correction_for_momentum = correction_for_heat
    elif richardson_number < -0.01:
        # ----------------------
        # unstable
        # ----------------------
        x = (1.0 - 16.0 * (measurement_height - zero_displacement_height) / monin_obukhov_length) ** 0.25
        y = x ** 2
        correction_for_heat = 2.0 * log((1 + y) / 2)  # (Liu et al. 2007, eq. 13)
        # Colaizzi et al. 2004, eq. 12)
        correction_for_momentum = 2.0 * log((1 + x) / 2) + log((1 + x ** 2) / 2) - 2 * atan(x) + 2 * atan(1)
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

    return correction_for_momentum, correction_for_heat, richardson_number, monin_obukhov_length
