from crop_energy_balance.formalisms import lumped_leaves
from crop_energy_balance.utils import is_almost_equal, assert_trend


def test_calc_leaf_layer_boundary_conductance_to_vapor():
    assert 0 == lumped_leaves.calc_leaf_layer_boundary_conductance_to_vapor(wind_speed_at_canopy_height=0,
                                                                            upper_cumulative_leaf_area_index=0,
                                                                            lower_cumulative_leaf_area_index=1,
                                                                            wind_speed_extinction_coefficient=0.5,
                                                                            characteristic_length=0.01,
                                                                            shape_parameter=0.01)

    assert is_almost_equal(
        actual=lumped_leaves.calc_leaf_layer_boundary_conductance_to_vapor(wind_speed_at_canopy_height=1,
                                                                           upper_cumulative_leaf_area_index=0,
                                                                           lower_cumulative_leaf_area_index=1.e9,
                                                                           wind_speed_extinction_coefficient=2,
                                                                           characteristic_length=1,
                                                                           shape_parameter=1 / 3600.),
        desired=1, decimal=3)


def test_calc_leaf_layer_boundary_resistance_to_heat():
    try:
        lumped_leaves.calc_leaf_layer_boundary_resistance_to_heat(wind_speed_at_canopy_height=0,
                                                                  upper_cumulative_leaf_area_index=0,
                                                                  lower_cumulative_leaf_area_index=1,
                                                                  wind_speed_extinction_coefficient=0.5,
                                                                  characteristic_length=0.01,
                                                                  shape_parameter=0.01)
    except ZeroDivisionError as e:
        pass

    assert is_almost_equal(
        actual=lumped_leaves.calc_leaf_layer_boundary_resistance_to_heat(wind_speed_at_canopy_height=1,
                                                                         upper_cumulative_leaf_area_index=0,
                                                                         lower_cumulative_leaf_area_index=1.e9,
                                                                         wind_speed_extinction_coefficient=2,
                                                                         characteristic_length=1,
                                                                         shape_parameter=1 / 3600.),
        desired=1, decimal=3)


def test_calc_leaf_layer_surface_conductance_to_vapor():
    def set_args():
        return dict(incident_direct_irradiance=100,
                    incident_diffuse_irradiance=100,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    stomatal_sensibility_to_water_status=1,
                    leaf_scattering_coefficient=0.15,
                    canopy_reflectance_to_direct_irradiance=0.5,
                    canopy_reflectance_to_diffuse_irradiance=0.5,
                    direct_extinction_coefficient=1,
                    direct_black_extinction_coefficient=1,
                    diffuse_extinction_coefficient=1,
                    maximum_stomatal_conductance=40,
                    residual_stomatal_conductance=0.4,
                    shape_parameter=105,
                    sublayers_number=5)

    args = set_args()
    args.update({'incident_direct_irradiance': 0, 'incident_diffuse_irradiance': 0})
    assert is_almost_equal(actual=lumped_leaves.calc_leaf_layer_surface_conductance_to_vapor(**args), desired=0.4)

    args = set_args()
    g = []
    for par in range(0, 300, 50):
        args.update({'incident_direct_irradiance': par, 'incident_diffuse_irradiance': par})
        g.append(lumped_leaves.calc_leaf_layer_surface_conductance_to_vapor(**args))
    assert_trend(values=g, expected_trend='+')

    args = set_args()
    g = []
    for lai in range(10):
        args.update({'lower_cumulative_leaf_area_index': lai})
        g.append(lumped_leaves.calc_leaf_layer_surface_conductance_to_vapor(**args))
    assert_trend(values=g, expected_trend='+')


def test_calc_leaf_layer_surface_resistance_to_vapor():
    def set_args():
        return dict(incident_direct_irradiance=100,
                    incident_diffuse_irradiance=100,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    stomatal_sensibility_to_water_status=1,
                    leaf_scattering_coefficient=0.15,
                    canopy_reflectance_to_direct_irradiance=0.5,
                    canopy_reflectance_to_diffuse_irradiance=0.5,
                    direct_extinction_coefficient=1,
                    direct_black_extinction_coefficient=1,
                    diffuse_extinction_coefficient=1,
                    maximum_stomatal_conductance=40,
                    residual_stomatal_conductance=0.4,
                    shape_parameter=105,
                    sublayers_number=5,
                    stomatal_density_factor=1)

    args = set_args()
    args.update({'residual_stomatal_conductance': 0, 'incident_direct_irradiance': 0, 'incident_diffuse_irradiance': 0})
    assert 1.e6 == lumped_leaves.calc_leaf_layer_surface_resistance_to_vapor(**args)

    args = set_args()
    args.update({'residual_stomatal_conductance': 1.e9})
    assert is_almost_equal(actual=lumped_leaves.calc_leaf_layer_surface_resistance_to_vapor(**args), desired=0)


def test_calc_leaf_layer_net_longwave_radiation():
    assert 0 == lumped_leaves.calc_leaf_layer_net_longwave_radiation(canopy_top_net_longwave_radiation=0,
                                                                     upper_cumulative_leaf_area_index=0,
                                                                     lower_cumulative_leaf_area_index=1,
                                                                     diffuse_black_extinction_coefficient=0.5)

    assert_trend(values=[lumped_leaves.calc_leaf_layer_net_longwave_radiation(canopy_top_net_longwave_radiation=1,
                                                                              upper_cumulative_leaf_area_index=0,
                                                                              lower_cumulative_leaf_area_index=lai,
                                                                              diffuse_black_extinction_coefficient=0.5)
                         for lai in range(10)],
                 expected_trend='+')
