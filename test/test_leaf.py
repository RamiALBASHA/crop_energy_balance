from crop_energy_balance.formalisms import leaf
from crop_energy_balance.utils import is_almost_equal, assert_trend


def test_calc_leaf_boundary_conductance():
    assert 0 == leaf.calc_forced_convection_conductance(wind_speed_at_canopy_height=0,
                                                        characteristic_length=1,
                                                        shape_parameter=0.01)

    assert is_almost_equal(actual=leaf.calc_forced_convection_conductance(wind_speed_at_canopy_height=1,
                                                                          characteristic_length=1.e9,
                                                                          shape_parameter=0.01),
                           desired=0, decimal=3)

    assert 36 == leaf.calc_forced_convection_conductance(wind_speed_at_canopy_height=1,
                                                         characteristic_length=1,
                                                         shape_parameter=0.01)


def test_calc_stomatal_sensibility_leuning():
    assert 1 == leaf.calc_stomatal_sensibility_leuning(vapor_pressure_deficit=0, d_0=1)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_leuning(vapor_pressure_deficit=1, d_0=1.e9),
        desired=1,
        decimal=3)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_leuning(vapor_pressure_deficit=1.e9, d_0=1),
        desired=0,
        decimal=3)

    assert_trend(
        values=[leaf.calc_stomatal_sensibility_leuning(vpd, 1) for vpd in range(7)],
        expected_trend='-')

    pass


def test_calc_stomatal_sensibility_tuzet():
    assert 1 == leaf.calc_stomatal_sensibility_tuzet(water_potential=0, psi_ref=0, steepness=0)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_tuzet(water_potential=0, psi_ref=0, steepness=-1.e9),
        desired=1,
        decimal=3)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_tuzet(water_potential=-1.e3, psi_ref=-300, steepness=1),
        desired=0,
        decimal=3)

    assert_trend(
        values=[leaf.calc_stomatal_sensibility_tuzet(-psi, -300, 1) for psi in range(0, 1000, 100)],
        expected_trend='-')

    pass


def test_calc_stomatal_sensibility_misson():
    assert 1 == leaf.calc_stomatal_sensibility_misson(water_potential=0, psi_half_aperture=-300, steepness=1)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_misson(water_potential=-600, psi_half_aperture=-300, steepness=-1.e3),
        desired=1.0,
        decimal=3)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_misson(water_potential=-300, psi_half_aperture=-300, steepness=-1.e9),
        desired=0.5,
        decimal=3)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility_misson(water_potential=-1.e6, psi_half_aperture=-300, steepness=1),
        desired=0,
        decimal=3)

    assert_trend(
        values=[leaf.calc_stomatal_sensibility_misson(psi, -300, 2) for psi in range(1, 1000, 100)],
        expected_trend='-')

    pass


def test_calc_stomatal_conductance():
    g_res = 4
    g_max = 40
    assert g_res == leaf.calc_stomatal_conductance(residual_stomatal_conductance=g_res,
                                                   maximum_stomatal_conductance=g_max,
                                                   absorbed_irradiance=0,
                                                   shape_parameter=1,
                                                   stomatal_sensibility_to_water_status=1)
    assert g_res + g_max / 2. == leaf.calc_stomatal_conductance(residual_stomatal_conductance=g_res,
                                                                maximum_stomatal_conductance=g_max,
                                                                absorbed_irradiance=1,
                                                                shape_parameter=1,
                                                                stomatal_sensibility_to_water_status=1)

    assert is_almost_equal(actual=leaf.calc_stomatal_conductance(residual_stomatal_conductance=g_res,
                                                                 maximum_stomatal_conductance=g_max,
                                                                 absorbed_irradiance=1,
                                                                 shape_parameter=1.e9,
                                                                 stomatal_sensibility_to_water_status=1),
                           desired=g_res, decimal=3)
