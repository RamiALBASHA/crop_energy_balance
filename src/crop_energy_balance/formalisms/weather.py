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


def calc_atmospheric_emissivity(air_vapor_pressure: float,
                                air_temperature: float) -> float:
    """Calculates the effective atmospheric emissivity of a clear sky

    Args:
        air_vapor_pressure: [kPa] air vapor pressure
        air_temperature: [K] air temperature

    Returns:
        [-] effective atmospheric emissivity of a clear sky

    References:
        Brutsaert, 1975.
            On a Derivable Formula for Long-Wave Radiation From Clear Skies.
            Water Resources Research 11, 742 - 744.
    """

    return 1.24 * (0.1 * air_vapor_pressure / air_temperature) ** (1. / 7.)


def calc_vapor_pressure_slope(temperature: float) -> float:
    """Calculates the slope of air vapor pressure curve at a given air temperature.

    Args:
        temperature: [°C] air temperature

    Returns:
        [kPa °K-1] the slope of vapor pressure curve at the given air temperature

    """
    return 4098 * (0.6108 * exp((17.27 * temperature) / (temperature + 237.3))) / ((temperature + 237.3) ** 2)
