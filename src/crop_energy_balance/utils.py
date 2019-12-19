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


def discretize_linearly(inclusive_start: float,
                        inclusive_stop: float,
                        vector_length: int) -> list:
    """Discretizes linearily an axis.

    Args:
        inclusive_start: inclusive start
        inclusive_stop: inclusive stop
        vector_length: desired length of the returned list

    Returns:
        Discretized axis

    Notes:
        This function replaces numpy.linspace()
    """
    assert vector_length != 1, 'args:`vector_length` must be greater than 1.'

    step = (inclusive_stop - inclusive_start) / (vector_length - 1)
    return [inclusive_start + step * i for i in range(vector_length)]
