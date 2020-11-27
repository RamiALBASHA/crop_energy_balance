from math import exp


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
        return 1
    else:
        return 2


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
        This function is intended to replaces numpy.linspace()
    """
    assert vector_length != 1, 'args:`vector_length` must be greater than 1.'

    step = (inclusive_stop - inclusive_start) / (vector_length - 1)
    return [inclusive_start + step * i for i in range(vector_length)]


def calc_temperature_step(previous_value: float,
                          actual_value: float,
                          step_fraction: float = 0.5) -> float:
    """Calculates the temperature step value between two consecutive energy balance calculations.

    Args:
        previous_value: [K] previous calculated temperature
        actual_value: [K] actual calculated temperature
        step_fraction: [-] fraction of the entire step (`actual_value - previous_value`) to be used

    Returns:
        [K] the temperature step value between two consecutive energy balance calculations
    """
    return step_fraction * (actual_value - previous_value)


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
