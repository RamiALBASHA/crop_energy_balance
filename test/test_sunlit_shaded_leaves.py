from crop_energy_balance.formalisms import sunlit_shaded_leaves
from crop_energy_balance.utils import is_almost_equal, assert_trend


def test_calc_leaf_layer_boundary_conductance_to_vapor():
    def set_args(**kwargs):
        args = dict(wind_speed_at_canopy_height=3600,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    direct_black_extinction_coefficient=0.5,
                    wind_speed_extinction_coefficient=0.5,
                    characteristic_length=0.01,
                    shape_parameter=0.01)
        args.update(**kwargs)
        return args

    for category, res in (('sunlit', 0), ('shaded', 1.e-12)):
        assert res == sunlit_shaded_leaves.calc_leaf_layer_boundary_conductance_to_vapor(
            **set_args(leaves_category=category, wind_speed_at_canopy_height=0))

        assert is_almost_equal(desired=0, actual=sunlit_shaded_leaves.calc_leaf_layer_boundary_conductance_to_vapor(
            **set_args(leaves_category=category,
                       upper_cumulative_leaf_area_index=1, lower_cumulative_leaf_area_index=1)))

        assert_trend(expected_trend='-',
                     values=[sunlit_shaded_leaves.calc_leaf_layer_boundary_conductance_to_vapor(
                         **set_args(leaves_category=category, wind_speed_extinction_coefficient=k_u))
                         for k_u in (0.1, 10, 100)])


def test_calc_leaf_layer_boundary_resistance_to_heat():
    def set_args(**kwargs):
        args = dict(wind_speed_at_canopy_height=3600,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    direct_black_extinction_coefficient=0.5,
                    wind_speed_extinction_coefficient=0.5,
                    characteristic_length=0.01,
                    shape_parameter=0.01)
        args.update(**kwargs)
        return args

    for category in ('sunlit', 'shaded'):
        try:
            sunlit_shaded_leaves.calc_leaf_layer_boundary_resistance_to_heat(
                **set_args(leaves_category=category, wind_speed_at_canopy_height=0))
        except ZeroDivisionError:
            pass

        assert_trend(expected_trend='+',
                     values=[sunlit_shaded_leaves.calc_leaf_layer_boundary_resistance_to_heat(
                         **set_args(leaves_category=category, wind_speed_extinction_coefficient=k_u))
                         for k_u in (0.1, 10, 100)])


def test_calc_leaf_layer_surface_conductance_to_vapor():
    def set_args(**kwargs):
        args = dict(incident_direct_irradiance=500,
                    incident_diffuse_irradiance=50,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    stomatal_sensibility_to_water_status=1,
                    leaf_scattering_coefficient=0.15,
                    canopy_reflectance_to_direct_irradiance=0.028,
                    canopy_reflectance_to_diffuse_irradiance=0.057,
                    direct_extinction_coefficient=0.5,
                    direct_black_extinction_coefficient=0.54,
                    diffuse_extinction_coefficient=0.62,
                    maximum_stomatal_conductance=39.6,
                    residual_stomatal_conductance=4,
                    shape_parameter=105,
                    sublayers_number=5)
        args.update(**kwargs)
        return args

    for category in ('sunlit', 'shaded'):
        assert 0 == sunlit_shaded_leaves.calc_leaf_layer_surface_conductance_to_vapor(
            **set_args(leaves_category=category,
                       residual_stomatal_conductance=0, incident_direct_irradiance=0, incident_diffuse_irradiance=0))

        assert_trend(expected_trend='+',
                     values=[sunlit_shaded_leaves.calc_leaf_layer_surface_conductance_to_vapor(**set_args(
                         leaves_category=category, incident_direct_irradiance=f, incident_diffuse_irradiance=f))
                             for f in range(1, 301, 50)])

        assert_trend(expected_trend='+',
                     values=[sunlit_shaded_leaves.calc_leaf_layer_surface_conductance_to_vapor(
                         **set_args(leaves_category=category, lower_cumulative_leaf_area_index=lai))
                         for lai in range(10)])


def test_calc_leaf_layer_surface_resistance_to_vapor():
    def set_args(**kwargs):
        args = dict(incident_direct_irradiance=500,
                    incident_diffuse_irradiance=50,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    stomatal_sensibility_to_water_status=1,
                    leaf_scattering_coefficient=0.15,
                    canopy_reflectance_to_direct_irradiance=0.028,
                    canopy_reflectance_to_diffuse_irradiance=0.057,
                    direct_extinction_coefficient=0.5,
                    direct_black_extinction_coefficient=0.54,
                    diffuse_extinction_coefficient=0.62,
                    maximum_stomatal_conductance=39.6,
                    residual_stomatal_conductance=4,
                    shape_parameter=105,
                    sublayers_number=5,
                    stomatal_density_factor=1)
        args.update(**kwargs)
        return args

    for category in ('sunlit', 'shaded'):
        try:
            sunlit_shaded_leaves.calc_leaf_layer_surface_resistance_to_vapor(
                **set_args(leaves_category=category, residual_stomatal_conductance=0,
                           incident_direct_irradiance=0, incident_diffuse_irradiance=0))
        except ZeroDivisionError:
            pass

        assert_trend(expected_trend='-',
                     values=[sunlit_shaded_leaves.calc_leaf_layer_surface_resistance_to_vapor(**set_args(
                         leaves_category=category, incident_direct_irradiance=f, incident_diffuse_irradiance=f))
                             for f in range(1, 301, 50)])

        assert_trend(expected_trend='-',
                     values=[sunlit_shaded_leaves.calc_leaf_layer_surface_resistance_to_vapor(
                         **set_args(leaves_category=category, lower_cumulative_leaf_area_index=lai))
                         for lai in range(10)])


def test_calc_leaf_fraction_per_leaf_layer():
    def set_args(**kwargs):
        args = dict(upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1.e-6,
                    direct_black_extinction_coefficient=0.5)
        args.update(**kwargs)
        return args

    assert is_almost_equal(
        desired=1,
        actual=(sunlit_shaded_leaves.calc_leaf_fraction_per_leaf_layer(**set_args(leaves_category='sunlit')) +
                sunlit_shaded_leaves.calc_leaf_fraction_per_leaf_layer(**set_args(leaves_category='shaded'))))


def test_calc_leaf_layer_net_longwave_radiation():
    def set_args(**kwargs):
        args = dict(canopy_top_net_longwave_radiation=-50,
                    upper_cumulative_leaf_area_index=0,
                    lower_cumulative_leaf_area_index=1,
                    direct_black_extinction_coefficient=0.54,
                    diffuse_black_extinction_coefficient=0.6)
        args.update(**kwargs)
        return args

    category = 'sunlit'

    assert 0 == sunlit_shaded_leaves.calc_leaf_layer_net_longwave_radiation(
        **set_args(leaves_category=category, canopy_top_net_longwave_radiation=0))

    assert_trend(expected_trend='-',
                 values=[sunlit_shaded_leaves.calc_leaf_layer_net_longwave_radiation(**set_args(
                     leaves_category=category, lower_cumulative_leaf_area_index=lai)) for lai in range(0, 10)])
