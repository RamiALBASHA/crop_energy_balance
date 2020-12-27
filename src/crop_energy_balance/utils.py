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
        This function is intended to replace numpy.linspace()
    """
    assert vector_length != 1, 'args:`vector_length` must be greater than 1.'

    step = (inclusive_stop - inclusive_start) / (vector_length - 1)
    return [inclusive_start + step * i for i in range(vector_length)]


def assert_almost_equal(actual: any([int, float, tuple, list]),
                        desired: any([int, float, tuple, list]),
                        decimal: int = 7):
    """Raises an AssertionError if two items are not equal up to desired precision.

    Args:
        actual: The object to check
        desired: The expected object
        decimal: Desired precision

    Raises:
        AssertionError if actual and desired are not equal up to specified precision.

    Notes:
        This function is a simplified version of numpy.testing.assert_almost_equal() and intends to replace the latter.
        The test verifies that the elements of ``actual`` and ``desired`` satisfy.
            ``abs(desired-actual) < 1.5 * 10**(-decimal)``

    """
    __tracebackhide__ = True  # Hide traceback for py.test

    try:
        any([isinstance(item, (tuple, list)) for item in (actual, desired)])
        assert len(actual) == len(desired)
        return [assert_almost_equal(actual=a, desired=d) for a, d in zip(actual, desired)]
    except TypeError:
        pass

    if abs(desired - actual) >= 1.5 * 10.0 ** (-decimal):
        raise AssertionError()


def assert_trend(values: list, expected_trend: str) -> None or AssertionError:
    """Asserts that a vector of values follows a given trend.

    Args:
        values: values whose trend is to be checked
        expected_trend: one of '+' (increasing), '-' (decreasing), '+-' (non-monotonic)
    """
    if expected_trend == '+':
        assert all([x <= y for x, y in zip(values, values[1:])])
    elif expected_trend == '-':
        assert all([x >= y for x, y in zip(values, values[1:])])
    elif expected_trend == '+-':
        assert (not all([x <= y for x, y in zip(values, values[1:])]) and
                not all([x >= y for x, y in zip(values, values[1:])]))
