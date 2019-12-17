def kelvin_to_celsius(temperature, absolute_zero=-273.15):
    """Converts Kelvin unit to Celsius unit.

    Args:
        temperature (float): [K] temperature
        absolute_zero (float): [°C] absolute zero temperature

    Returns:
        (float): [K] temperature

    """
    return temperature + absolute_zero


def celsius_to_kelvin(temperature, absolute_zero=-273.15):
    """Converts Kelvin unit to Celsius unit.

    Args:
        temperature (float): [°C] temperature
        absolute_zero (float): [°C] absolute zero temperature

    Returns:
        (float): [°C] temperature

    """
    return temperature - absolute_zero


def stomatal_density_factor(amphistomatal_leaf):
    """Computes an integer that expresses whether stomata are equally present on both faces or on one face of the leaf
        blade

    Args:
        amphistomatal_leaf (bool): `True` if both blade sides have stomata of equivalent density

    Returns:
        (int): [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    """
    if amphistomatal_leaf:
        nu = 1
    else:
        nu = 2

    return nu
