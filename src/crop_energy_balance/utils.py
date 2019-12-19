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


def calc_stomatal_density_factor(amphistomatal_leaf: bool) -> int:
    """Computes an integer that expresses whether stomata are equally present on both faces or on one face of the leaf
        blade

    Args:
        amphistomatal_leaf: `True` if both blade sides have stomata of equivalent density

    Returns:
        [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    """
    if amphistomatal_leaf:
        nu = 1
    else:
        nu = 2

    return nu
