from math import exp


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
        [kPa °K-1] the slope of vapor pressure curve at the given air temperature

    """
    return 4098 * (0.6108 * exp((17.27 * temperature) / (temperature + 237.3))) / ((temperature + 237.3) ** 2)


def convert_kelvin_to_celsius(temperature: float,
                              absolute_zero: float = -273.15) -> float:
    """Converts Kelvin unit to Celsius unit.

    Args:
        temperature: [K] temperature
        absolute_zero: [°C] absolute zero temperature

    Returns:
        [K] temperature

    """
    return temperature + absolute_zero


def convert_celsius_to_kelvin(temperature: float,
                              absolute_zero: float = -273.15) -> float:
    """Converts Kelvin unit to Celsius unit.

    Args:
        temperature: [°C] temperature
        absolute_zero: [°C] absolute zero temperature

    Returns:
        [°C] temperature

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
        relative_humidity: [-] air relative humidity (%, between 0 and 1)
    Returns:
        [kPa] leaf-to-air vapour pressure deficit
    """
    es_l = calc_saturated_air_vapor_pressure(temperature_leaf)
    es_a = calc_saturated_air_vapor_pressure(temperature_air)
    ea = es_a * relative_humidity / 100

    return es_l - ea
