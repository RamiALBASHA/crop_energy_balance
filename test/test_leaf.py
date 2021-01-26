from crop_energy_balance.formalisms import leaf
from crop_energy_balance.utils import is_almost_equal


def test_calc_leaf_boundary_conductance():
    assert 0 == leaf.calc_leaf_boundary_conductance(wind_speed_at_canopy_height=0,
                                                    characteristic_length=1,
                                                    shape_parameter=0.01)

    assert is_almost_equal(actual=leaf.calc_leaf_boundary_conductance(wind_speed_at_canopy_height=1,
                                                                      characteristic_length=1.e9,
                                                                      shape_parameter=0.01),
                           desired=0, decimal=3)

    assert 36 == leaf.calc_leaf_boundary_conductance(wind_speed_at_canopy_height=1,
                                                     characteristic_length=1,
                                                     shape_parameter=0.01)


def test_calc_stomatal_sensibility():
    assert 1 == leaf.calc_stomatal_sensibility(air_vapor_pressure_deficit=0, shape_parameter=1)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility(air_vapor_pressure_deficit=1, shape_parameter=1.e9), desired=1, decimal=3)

    assert is_almost_equal(
        actual=leaf.calc_stomatal_sensibility(air_vapor_pressure_deficit=1.e9, shape_parameter=1), desired=0, decimal=3)

    try:
        leaf.calc_stomatal_sensibility(air_vapor_pressure_deficit=0, shape_parameter=0)
    except AssertionError as e:
        assert e.args[0] == 'The value of `shape_parameter` must be greater than zero.'


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
