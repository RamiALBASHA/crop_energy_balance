from crop_energy_balance.formalisms import soil
from crop_energy_balance.params import Constants
from crop_energy_balance.utils import is_almost_equal, assert_trend

constants = Constants()


def test_calc_eddy_diffusivity():
    u = 2400
    z_h = 1
    z_r = 2
    d = 0.67
    z0 = 0.123
    k = constants.von_karman

    assert 0 == soil.calc_eddy_diffusivity(wind_speed=0,
                                           canopy_height=z_h,
                                           measurement_height=z_r,
                                           zero_displacement_height=d,
                                           canopy_roughness_length_for_momentum=z0,
                                           von_karman_constant=k)

    assert 0 == soil.calc_eddy_diffusivity(wind_speed=u,
                                           canopy_height=d,
                                           measurement_height=z_r,
                                           zero_displacement_height=d,
                                           canopy_roughness_length_for_momentum=z0,
                                           von_karman_constant=k)

    assert is_almost_equal(desired=1,
                           actual=soil.calc_eddy_diffusivity(wind_speed=1 / k ** 2,
                                                             canopy_height=d + 1,
                                                             measurement_height=2 * d,
                                                             zero_displacement_height=d,
                                                             canopy_roughness_length_for_momentum=d / 2.718281828459045,
                                                             von_karman_constant=k))

    try:
        soil.calc_eddy_diffusivity(wind_speed=1 / k ** 2,
                                   canopy_height=d + 1,
                                   measurement_height=2 * d,
                                   zero_displacement_height=d,
                                   canopy_roughness_length_for_momentum=d,
                                   von_karman_constant=k)
    except ZeroDivisionError:
        pass


def test_calc_boundary_resistance():
    def set_args(**kwargs):
        args = dict(wind_speed=2400,
                    canopy_height=1,
                    measurement_height=2,
                    zero_displacement_height=0.67,
                    canopy_roughness_length_for_momentum=0.123,
                    soil_roughness_length_for_momentum=0.01,
                    shape_parameter=2.5,
                    von_karman_constant=constants.von_karman)
        args.update(**kwargs)
        return args

    args_test_1 = set_args(soil_roughness_length_for_momentum=1, canopy_roughness_length_for_momentum=1,
                           zero_displacement_height=0, shape_parameter=1)
    assert 0 == soil.calc_forced_convection_resistance(**args_test_1)

    args_test_2 = set_args(wind_speed=0)
    try:
        soil.calc_forced_convection_resistance(**args_test_2)
    except ZeroDivisionError:
        pass


def test_calc_surface_resistance():
    assert 1. / 3600 == soil.calc_surface_resistance(soil_saturation_ratio=1.,
                                                     shape_parameter_1=1.,
                                                     shape_parameter_2=1.)

    assert_trend(expected_trend='+',
                 values=[soil.calc_surface_resistance(soil_saturation_ratio=sat,
                                                      shape_parameter_1=1.,
                                                      shape_parameter_2=1.)
                         for sat in (1, 0.5, 0)])


def test_calc_heat_flux():
    rn = 100

    assert 10 == soil.calc_heat_flux(net_above_ground_radiation=rn, is_diurnal=True)

    assert 50 == soil.calc_heat_flux(net_above_ground_radiation=rn, is_diurnal=False)


def test_calc_net_longwave_radiation():
    rnl = 50
    lai = 3
    k_d = 1

    assert 0 == soil.calc_net_longwave_radiation(canopy_top_net_longwave_radiation=0,
                                                 canopy_leaf_area_index=lai,
                                                 diffuse_black_extinction_coefficient=k_d)

    assert rnl == soil.calc_net_longwave_radiation(canopy_top_net_longwave_radiation=rnl,
                                                   canopy_leaf_area_index=0,
                                                   diffuse_black_extinction_coefficient=k_d)

    assert_trend(expected_trend='-',
                 values=[soil.calc_net_longwave_radiation(canopy_top_net_longwave_radiation=rnl,
                                                          canopy_leaf_area_index=lai,
                                                          diffuse_black_extinction_coefficient=k_d)
                         for lai in range(0, 5)])
