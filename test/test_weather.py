from crop_energy_balance.formalisms import weather
from crop_energy_balance.params import Constants
from crop_energy_balance.utils import is_almost_equal, assert_trend

constants = Constants()


def test_calc_psychrometric_constant():
    assert_trend(values=[weather.calc_psychrometric_constant(atmospheric_pressure=p) for p in range(90, 102, 2)],
                 expected_trend='+')


def test_calc_atmospheric_emissivity():
    p_air = 3
    t_air = 25 + 273
    try:
        weather.calc_atmospheric_emissivity(model='some unknown model', air_vapor_pressure=3, air_temperature=25 + 273)
    except ValueError as e:
        assert e.args[0].startswith('Unknown model name')

    for model in ('brutsaert_1975', 'monteith_2013'):
        assert_trend(
            expected_trend='+',
            values=[weather.calc_atmospheric_emissivity(model=model, air_vapor_pressure=p, air_temperature=t_air)
                    for p in range(0, 7)])

        assert_trend(
            expected_trend='-',
            values=[weather.calc_atmospheric_emissivity(model=model, air_vapor_pressure=p_air, air_temperature=t)
                    for t in range(-25 + 273, 25 + 273, 10)])


def test_calc_vapor_pressure_slope():
    assert_trend(expected_trend='+',
                 values=[weather.calc_vapor_pressure_slope(t) for t in range(-25, 25)])
    # (Allen et al., 1998, Annex 2, table 2.4)
    air_temperature = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45]
    slope = [0.047, 0.061, 0.082, 0.11, 0.145, 0.189, 0.243, 0.311, 0.393, 0.493]
    assert is_almost_equal(
        actual=[weather.calc_vapor_pressure_slope(t) for t in air_temperature], desired=slope, decimal=3)


def test_convert_kelvin_to_celsius():
    assert constants.absolute_zero == weather.convert_kelvin_to_celsius(temperature=0,
                                                                        absolute_zero=constants.absolute_zero)


def test_convert_celsius_to_kelvin():
    assert -constants.absolute_zero == weather.convert_celsius_to_kelvin(temperature=0,
                                                                         absolute_zero=constants.absolute_zero)


def test_convert_photosynthetically_active_radiation_into_global_radiation():
    assert 1 / 0.48 == weather.convert_photosynthetically_active_radiation_into_global_radiation(1)


def test_convert_global_irradiance_into_photosynthetically_active_radiation():
    assert 0.48 == weather.convert_global_irradiance_into_photosynthetically_active_radiation(1)


def test_calc_saturated_air_vapor_pressure():
    assert_trend(expected_trend='+',
                 values=[weather.calc_saturated_air_vapor_pressure(t) for t in range(-25, 25)])


def test_calc_vapor_pressure_deficit():
    assert 0 == weather.calc_vapor_pressure_deficit(temperature_air=25, temperature_leaf=25, relative_humidity=100)

    assert_trend(
        expected_trend='-',
        values=[weather.calc_vapor_pressure_deficit(temperature_air=t, temperature_leaf=25, relative_humidity=100)
                for t in range(25, 50)])

    assert_trend(
        expected_trend='+',
        values=[weather.calc_vapor_pressure_deficit(temperature_air=25, temperature_leaf=t, relative_humidity=100)
                for t in range(25, 50)])

    assert_trend(
        expected_trend='-',
        values=[weather.calc_vapor_pressure_deficit(temperature_air=25, temperature_leaf=25, relative_humidity=rh)
                for rh in range(0, 100)])
