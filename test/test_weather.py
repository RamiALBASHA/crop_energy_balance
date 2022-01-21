import crop_energy_balance.formalisms.canopy
from crop_energy_balance.formalisms import weather
from crop_energy_balance.params import Constants
from crop_energy_balance.utils import is_almost_equal, assert_trend, discretize_linearly

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


def test_calc_monin_obukhov_length():
    def set_args(**kwargs):
        args = dict(surface_temperature=25 + 273,
                    sensible_heat_flux=50,
                    friction_velocity=400,
                    air_density=constants.air_density,
                    air_specific_heat_capacity=constants.air_specific_heat_capacity,
                    gravitational_acceleration=constants.gravitational_acceleration,
                    von_karman_constant=constants.von_karman)
        args.update(**kwargs)
        return args

    assert 0 == weather.calc_monin_obukhov_length(**set_args(friction_velocity=0))

    assert_trend(expected_trend='-',
                 values=[weather.calc_monin_obukhov_length(**set_args(surface_temperature=t)) for t in range(248, 298)])

    assert_trend(expected_trend='+',
                 values=[weather.calc_monin_obukhov_length(**set_args(sensible_heat_flux=h)) for h in range(1, 10, 90)])

    assert_trend(expected_trend='-',
                 values=[weather.calc_monin_obukhov_length(**set_args(sensible_heat_flux=h)) for h in range(-1, -10)])


def test_calc_richardson_number():
    def set_args(**kwargs):
        args = dict(measurement_height=2, zero_displacement_height=0.67, monin_obukhov_length=5)
        args.update(**kwargs)
        return args

    assert (weather.calc_richardson_number(**set_args(is_stable=True)) !=
            weather.calc_richardson_number(**set_args(is_stable=False)))

    assert 0 == weather.calc_richardson_number(
        **set_args(is_stable=True, measurement_height=2, zero_displacement_height=2))

    assert 0 == weather.calc_richardson_number(
        **set_args(is_stable=False, measurement_height=2, zero_displacement_height=2))


def test_calc_stability_correction_functions():
    def set_args(**kwargs):
        args = dict(friction_velocity=400,
                    sensible_heat=100,
                    canopy_temperature=5 + 273,
                    measurement_height=2,
                    zero_displacement_height=0.67,
                    air_density=constants.air_density,
                    air_specific_heat_capacity=constants.air_specific_heat_capacity,
                    von_karman_constant=constants.von_karman,
                    gravitational_acceleration=constants.gravitational_acceleration)
        args.update(**kwargs)
        return args

    assert 4 == len(weather.calc_stability_correction_functions(**set_args()))

    phi_m, phi_h, richardson, obukhov = weather.calc_stability_correction_functions(**set_args(friction_velocity=0.1))
    assert is_almost_equal(actual=[phi_m, phi_h, obukhov], desired=[0, 0, 0])

    phi_m, phi_h, richardson, obukhov = weather.calc_stability_correction_functions(**set_args(friction_velocity=1.e6))
    assert is_almost_equal(actual=richardson, desired=0)
    assert phi_m != 0
    assert phi_h != 0
