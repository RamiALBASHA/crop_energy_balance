from crop_energy_balance import utils
from crop_energy_balance.formalisms import canopy, weather
from crop_energy_balance.params import Constants
from crop_energy_balance.utils import is_almost_equal, assert_trend, discretize_linearly

constants = Constants()


def test_calc_zero_displacement_height():
    assert canopy.calc_zero_displacement_height(0, 1, 0.2) == 0
    assert round(canopy.calc_zero_displacement_height(1, 1, 0.2), 3) == 0.563


def test_calc_roughness_length_for_momentum_transfer():
    canopy_height = 1
    drag_coefficient = 0.2
    values = [canopy.calc_roughness_length_for_momentum_transfer(
        soil_roughness_length_for_momentum=0.01,
        canopy_height=canopy_height,
        zero_plan_displacement_height=canopy.calc_zero_displacement_height(
            canopy_height=canopy_height,
            leaf_area_index=lai,
            drag_coefficient=drag_coefficient
        ),
        total_leaf_area_index=lai / 10,
        drag_coefficient=0.2) for lai in range(40)[1:]]
    assert_trend(values=values, expected_trend='+-')


def test_calc_roughness_length_for_heat_transfer():
    assert canopy.calc_roughness_length_for_heat_transfer(0, 0.1) == 0
    assert canopy.calc_roughness_length_for_heat_transfer(1, 0.1) == 0.1


def test_calc_wind_speed_at_canopy_height():
    assert (canopy.calc_wind_speed_at_canopy_height(wind_speed=2400, canopy_height=1, measurement_height=2,
                                                    zero_displacement_height=0.67,
                                                    roughness_length_for_momentum=0.123) <
            canopy.calc_wind_speed_at_canopy_height(wind_speed=2401, canopy_height=1, measurement_height=2,
                                                    zero_displacement_height=0.67,
                                                    roughness_length_for_momentum=0.123))
    # minimal canopy height set to 0.1 m
    assert (canopy.calc_wind_speed_at_canopy_height(wind_speed=3600, canopy_height=0, measurement_height=2,
                                                    zero_displacement_height=0.067,
                                                    roughness_length_for_momentum=0.0123) ==
            canopy.calc_wind_speed_at_canopy_height(wind_speed=3600, canopy_height=0.1, measurement_height=2,
                                                    zero_displacement_height=0.067,
                                                    roughness_length_for_momentum=0.0123))

    assert canopy.calc_wind_speed_at_canopy_height(wind_speed=3600, canopy_height=1, measurement_height=1,
                                                   zero_displacement_height=0.67,
                                                   roughness_length_for_momentum=0.123) == 3600

    assert canopy.calc_wind_speed_at_canopy_height(wind_speed=3600, canopy_height=1, measurement_height=2,
                                                   zero_displacement_height=0.67,
                                                   roughness_length_for_momentum=0.123) < 3600


def test_calc_net_longwave_radiation():
    vp = 3  # air_vapor_pressure
    t_air = weather.convert_celsius_to_kelvin(25)  # air_temperature
    atm = 0.73  # atmospheric_emissivity
    boltzmann = constants.stefan_boltzmann

    utils.assert_trend(
        values=[canopy.calc_net_longwave_radiation(weather.convert_celsius_to_kelvin(t), vp, atm, boltzmann)
                for t in range(-50, 50)],
        expected_trend='-')

    utils.assert_trend(
        values=[canopy.calc_net_longwave_radiation(t_air, vp, atm_emissivity / 100., boltzmann)
                for atm_emissivity in range(0, 110, 10)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_net_longwave_radiation(t_air, air_vp / 10., atm, boltzmann)
                for air_vp in range(0, 110, 10)],
        expected_trend='+')


def test_calc_sensible_heat_flux():
    t_source = weather.convert_celsius_to_kelvin(25)  # source_temperature
    t_air = weather.convert_celsius_to_kelvin(25)  # air_temperature
    r0 = 1  # aerodynamic_resistance
    rho = constants.air_density
    cp = constants.air_specific_heat_capacity

    assert canopy.calc_sensible_heat_flux(t_air, t_air, r0, rho, cp) == 0

    utils.assert_trend(
        values=[canopy.calc_sensible_heat_flux(weather.convert_celsius_to_kelvin(t), t_air, r0, rho, cp)
                for t in range(25, 50)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_sensible_heat_flux(t_source, t_air, r, rho, cp) for r in (0.1, 0.5, 1., 2., 3.)],
        expected_trend='+')


def test_calc_canopy_aerodynamic_resistance_under_neutral_conditions():
    """This test is taken from Box 4 in Allen et al. (1998). FAO Irrigation and Drainage Paper No. 56. Eq. 8
    """
    von_karman_constant = constants.von_karman
    height = 0.12
    measurement_height = 2
    leaf_area_index = 3
    drag_coef = 0.2
    zero_displacement_height = canopy.calc_zero_displacement_height(
        canopy_height=height, leaf_area_index=1, drag_coefficient=0.2)
    roughness_length_for_momentum = canopy.calc_roughness_length_for_momentum_transfer(
        soil_roughness_length_for_momentum=0.0125,
        zero_plan_displacement_height=zero_displacement_height,
        canopy_height=height,
        total_leaf_area_index=leaf_area_index,
        drag_coefficient=drag_coef)
    friction_velocity = canopy.calc_friction_velocity(
        wind_speed=3600,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        roughness_length_for_momentum=roughness_length_for_momentum,
        stability_correction_for_momentum=0,
        von_karman_constant=von_karman_constant)

    aerodynamic_resistance = canopy.calc_aerodynamic_resistance(
        richardson_number=1,
        friction_velocity=friction_velocity,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        roughness_length_for_heat=0.1 * roughness_length_for_momentum,
        stability_correction_for_heat=0,
        canopy_temperature=25,
        air_temperature=25,
        richardon_threshold_free_convection=-0.8,
        von_karman_constant=von_karman_constant,
        air_density=constants.air_density,
        air_specific_heat_capacity=constants.air_specific_heat_capacity,
        free_convection_shape_parameter=1.52)
    print(aerodynamic_resistance)
    assert utils.is_almost_equal(aerodynamic_resistance, 208 / 3600., decimal=2)


def test_calc_penman_evaporative_energy():
    r0 = 0.015  # canopy_aerodynamic_resistance
    a = 940  # canopy_available_energy
    s = 0.18  # vapor_pressure_slope
    vpd = 2.3  # vapor_pressure_deficit
    gamma = constants.psychrometric_constant
    rho = constants.air_density
    cp = constants.air_specific_heat_capacity

    utils.assert_trend(
        values=[canopy.calc_penman_evaporative_energy(r0_, a, s, vpd, gamma, rho, cp) for r0_ in (0.001, 0.1, 1)],
        expected_trend='-')

    utils.assert_trend(
        values=[canopy.calc_penman_evaporative_energy(r0, a_, s, vpd, gamma, rho, cp) for a_ in range(0, 1000, 100)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_penman_evaporative_energy(r0, a, s_, vpd, gamma, rho, cp) for s_ in range(0, 10)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_penman_evaporative_energy(r0, a, s, vpd_, gamma, rho, cp) for vpd_ in (0, 0.5, 1, 2, 3)],
        expected_trend='+')


def test_calc_lumped_aerodynamic_resistance():
    r0 = 0.015  # canopy_aerodynamic_resistance
    s = 0.18  # vapor_pressure_slope
    gamma = constants.psychrometric_constant

    utils.assert_trend(
        values=[canopy.calc_lumped_aerodynamic_resistance(r0_, s, gamma) for r0_ in (0.001, 0.1, 1)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_lumped_aerodynamic_resistance(r0, s_, gamma) for s_ in range(0, 10)],
        expected_trend='+')


def test_calc_penman_monteith_evaporative_energy():
    r0_lumped = 0.05  # canopy_lumped_aerodynamic_resistance
    penman = 920  # penman_evaporative_energy
    p = [2.2, 12.8]  # composed_boundary_and_surface_conductances
    a = [53.7, 790.9]  # available_energy_fluxes
    r = [0.026, 0.001]  # boundary_layer_resistances
    s = 0.184  # vapor_pressure_slope
    gamma = constants.psychrometric_constant

    utils.assert_trend(
        values=[canopy.calc_penman_monteith_evaporative_energy(r0_lumped_, penman, p, a, r, s, gamma)
                for r0_lumped_ in (0, 0.01, 0.05)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_penman_monteith_evaporative_energy(r0_lumped, penman_, p, a, r, s, gamma)
                for penman_ in range(0, 1000, 100)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_penman_monteith_evaporative_energy(r0_lumped, penman, p_, a, r, s, gamma)
                for p_ in ([[0, 0], [1, 1], [10, 10]])],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_penman_monteith_evaporative_energy(r0_lumped, penman, p, a_, r, s, gamma)
                for a_ in ([[0, 0], [10, 10], [100, 100], [1000, 1000]])],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_penman_monteith_evaporative_energy(r0_lumped, penman, p, a, r_, s, gamma)
                for r_ in ([[0, 0], [0.01, 0.01], [0.1, 0.1], [1, 1]])],
        expected_trend='+')


def test_calc_temperature():
    t_air = 297.65  # air_temperature
    r0 = 0.013  # canopy_aerodynamic_resistance
    a = 940  # canopy_available_energy
    pm = 736  # penman_monteith_evaporative_energy
    rho = constants.air_density
    cp = constants.air_specific_heat_capacity

    assert canopy.calc_temperature(t_air, r0, a, a, rho, cp) == t_air

    utils.assert_trend(
        values=[canopy.calc_temperature(weather.convert_celsius_to_kelvin(t_air_), r0, a, pm, rho, cp)
                for t_air_ in range(-50, 50)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_temperature(t_air, r0_, a, pm, rho, cp) for r0_ in (0.001, 0.1, 1)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_temperature(t_air, r0, a_, pm, rho, cp) for a_ in range(0, 1000, 100)],
        expected_trend='+')

    utils.assert_trend(
        values=[canopy.calc_temperature(t_air, r0, a, pm_, rho, cp) for pm_ in range(0, 1000, 100)],
        expected_trend='-')


def test_calc_friction_velocity():
    def set_args(**kwargs):
        args = dict(wind_speed=3600,
                    measurement_height=2,
                    zero_displacement_height=0.67,
                    roughness_length_for_momentum=0.123,
                    stability_correction_for_momentum=0,
                    von_karman_constant=constants.von_karman)
        args.update(**kwargs)
        return args

    assert 0 == canopy.calc_friction_velocity(**set_args(wind_speed=0))

    assert is_almost_equal(
        desired=-1,
        actual=canopy.calc_friction_velocity(**set_args(
            wind_speed=1 / constants.von_karman, zero_displacement_height=1, roughness_length_for_momentum=1,
            stability_correction_for_momentum=1)))

    assert_trend(expected_trend='+',
                 values=[canopy.calc_friction_velocity(**set_args(stability_correction_for_momentum=phi_m))
                         for phi_m in discretize_linearly(0.01, 2.38, 3)])  # 2.38075 = log((2 - 0.67) / 0.123)


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

    assert 0 == canopy.calc_monin_obukhov_length(**set_args(friction_velocity=0))

    assert_trend(expected_trend='-',
                 values=[canopy.calc_monin_obukhov_length(**set_args(surface_temperature=t)) for t in range(248, 298)])

    assert_trend(expected_trend='+',
                 values=[canopy.calc_monin_obukhov_length(**set_args(sensible_heat_flux=h)) for h in range(1, 10, 90)])

    assert_trend(expected_trend='-',
                 values=[canopy.calc_monin_obukhov_length(**set_args(sensible_heat_flux=h)) for h in range(-1, -10)])


def test_calc_richardson_number():
    def set_args(**kwargs):
        args = dict(measurement_height=2, zero_displacement_height=0.67, monin_obukhov_length=5)
        args.update(**kwargs)
        return args

    assert (canopy.calc_richardson_number(**set_args(is_stable=True)) !=
            canopy.calc_richardson_number(**set_args(is_stable=False)))

    assert 0 == canopy.calc_richardson_number(
        **set_args(is_stable=True, measurement_height=2, zero_displacement_height=2))

    assert 0 == canopy.calc_richardson_number(
        **set_args(is_stable=False, measurement_height=2, zero_displacement_height=2))


def test_calc_stability_correction_functions():
    friction_velocity = 400
    sensible_heat = 10
    canopy_temperature = 5 + 273
    measurement_height = 2
    zero_displacement_height = 0.67
    air_density = constants.air_density
    air_specific_heat_capacity = constants.air_specific_heat_capacity
    von_karman_constant = constants.von_karman
    gravitational_acceleration = constants.gravitational_acceleration

    monin_obukhov_length = canopy.calc_monin_obukhov_length(
        surface_temperature=canopy_temperature,
        sensible_heat_flux=sensible_heat,
        friction_velocity=friction_velocity,
        air_density=air_density,
        air_specific_heat_capacity=air_specific_heat_capacity,
        gravitational_acceleration=gravitational_acceleration,
        von_karman_constant=von_karman_constant)

    richardson_number = canopy.calc_richardson_number(
        is_stable=sensible_heat < 0,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        monin_obukhov_length=monin_obukhov_length)

    phi_m, phi_h = canopy.calc_stability_correction_functions(
        monin_obukhov_length=monin_obukhov_length,
        richardson_number=-1.1,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        richardon_threshold_free_convection=-0.8)
    assert phi_m == 0
    assert phi_h == 0

    phi_m, phi_h = canopy.calc_stability_correction_functions(
        monin_obukhov_length=monin_obukhov_length,
        richardson_number=0.2,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        richardon_threshold_free_convection=-0.8)
    assert phi_m == 0
    assert phi_h == 0

    phi_m, phi_h = canopy.calc_stability_correction_functions(
        monin_obukhov_length=monin_obukhov_length,
        richardson_number=richardson_number,
        measurement_height=measurement_height,
        zero_displacement_height=zero_displacement_height,
        richardon_threshold_free_convection=-0.8)
    assert phi_m != 0
    assert phi_h != 0
