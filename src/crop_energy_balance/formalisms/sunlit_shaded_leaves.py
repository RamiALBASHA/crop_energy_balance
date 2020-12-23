from math import exp

from crop_energy_balance.formalisms import lumped_leaves, leaf
from crop_energy_balance.formalisms.irradiance import calc_leaf_fraction, calc_absorbed_irradiance
from crop_energy_balance.utils import discretize_linearly
from crop_irradiance.uniform_crops.formalisms.sunlit_shaded_leaves import calc_sunlit_fraction_per_leaf_layer


def calc_leaf_layer_boundary_conductance_to_vapor(leaves_category: str,
                                                  wind_speed_at_canopy_height: float,
                                                  upper_cumulative_leaf_area_index: float,
                                                  lower_cumulative_leaf_area_index: float,
                                                  direct_black_extinction_coefficient: float,
                                                  wind_speed_extinction_coefficient: float = 0.5,
                                                  characteristic_length: float = 0.01,
                                                  shape_parameter: float = 0.01) -> float:
    """Calculates bulk layer boundary layer conductance for both sides of leaves blade.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        direct_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of direct (beam) irradiance
            through a canopy of black leaves
        wind_speed_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter

    Returns:
        [m h-1] bulk layer boundary layer conductance to water vapor for both sides of leaves blade
    """
    leaf_boundary_conductance = leaf.calc_leaf_boundary_conductance(wind_speed_at_canopy_height,
                                                                    characteristic_length,
                                                                    shape_parameter)
    lumped_extinction_coefficient = 0.5 * wind_speed_extinction_coefficient + direct_black_extinction_coefficient
    sunlit_layer_scaling_factor = 1.0 / lumped_extinction_coefficient * (
            exp(-lumped_extinction_coefficient * upper_cumulative_leaf_area_index) -
            exp(-lumped_extinction_coefficient * lower_cumulative_leaf_area_index))

    sunlit_layer_boundary_conductance = leaf_boundary_conductance * max(1.e-6, sunlit_layer_scaling_factor)

    if leaves_category == 'sunlit':
        return sunlit_layer_boundary_conductance
    elif leaves_category == 'shaded':
        lumped_layer_boundary_conductance = lumped_leaves.calc_leaf_layer_boundary_conductance_to_vapor(
            wind_speed_at_canopy_height,
            upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index,
            wind_speed_extinction_coefficient,
            characteristic_length,
            shape_parameter)
        return lumped_layer_boundary_conductance - sunlit_layer_boundary_conductance


def calc_leaf_layer_boundary_resistance_to_heat(leaves_category: str,
                                                wind_speed_at_canopy_height: float,
                                                upper_cumulative_leaf_area_index: float,
                                                lower_cumulative_leaf_area_index: float,
                                                direct_black_extinction_coefficient: float,
                                                wind_speed_extinction_coefficient: float,
                                                characteristic_length: float,
                                                shape_parameter: float) -> float:
    """Calculates the bulk leaf layer resistance to heat transfer.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        direct_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of direct (beam) irradiance
            through a canopy of black leaves
        wind_speed_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter

    Returns:
        [h m-1] bulk leaf layer resistance to heat transfer
    """

    return 1.0 / calc_leaf_layer_boundary_conductance_to_vapor(**locals())


def calc_leaf_layer_surface_conductance_to_vapor(leaves_category: str,
                                                 incident_direct_irradiance: float,
                                                 incident_diffuse_irradiance: float,
                                                 upper_cumulative_leaf_area_index: float,
                                                 lower_cumulative_leaf_area_index: float,
                                                 stomatal_sensibility_to_water_status: float,
                                                 leaf_scattering_coefficient: float,
                                                 canopy_reflectance_to_direct_irradiance: float,
                                                 canopy_reflectance_to_diffuse_irradiance: float,
                                                 direct_extinction_coefficient: float,
                                                 direct_black_extinction_coefficient: float,
                                                 diffuse_extinction_coefficient: float,
                                                 maximum_stomatal_conductance: float,
                                                 residual_stomatal_conductance: float,
                                                 shape_parameter: float = 105,
                                                 sublayers_number: int = 5) -> float:
    """Calculates the bulk surface conductance of a leaf layer for both sides of leaves blade.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        incident_direct_irradiance: [W m-2ground] incident direct photosynthetically active radiation above the canopy
        incident_diffuse_irradiance: [W m-2ground] incident diffuse irradiance at the top of the canopy
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index above the considered layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index below the considered layer
        stomatal_sensibility_to_water_status: [-] stomatal closure fraction due to water stress
        leaf_scattering_coefficient: [-] leaf scattering coefficient
        canopy_reflectance_to_direct_irradiance: [-] canopy reflectance to direct (beam) irradiance
        canopy_reflectance_to_diffuse_irradiance: [-] canopy reflectance to diffuse irradiance
        direct_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance
        direct_black_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam)
            irradiance for black leaves
        diffuse_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of diffuse irradiance
        maximum_stomatal_conductance: [m h-1] maximum stomatal conductance
        residual_stomatal_conductance: [m h-1] residual (minimum) stomatal conductance
        shape_parameter: [W m-2leaf] an empirical parameter to regulate the shape of stomatal conductance response
            to the absorbed photosynthetically active radiation (PAR)
        sublayers_number: number of sublayers that are used to performe numerical integral of the leaf-layer surface
            conductance equation

    Returns:
        [m h-1] bulk surface conductance of the leaf layer for both sides of leaves blade
    """
    sublayer_thickness = (lower_cumulative_leaf_area_index - upper_cumulative_leaf_area_index) / sublayers_number

    leaf_layer_surface_conductance = 0
    for cumulative_leaf_area_index in discretize_linearly(upper_cumulative_leaf_area_index,
                                                          lower_cumulative_leaf_area_index, sublayers_number):
        absorbed_irradiance = calc_absorbed_irradiance(leaves_category,
                                                       incident_direct_irradiance,
                                                       incident_diffuse_irradiance,
                                                       cumulative_leaf_area_index,
                                                       leaf_scattering_coefficient,
                                                       canopy_reflectance_to_direct_irradiance,
                                                       canopy_reflectance_to_diffuse_irradiance,
                                                       direct_extinction_coefficient,
                                                       direct_black_extinction_coefficient,
                                                       diffuse_extinction_coefficient)

        lumped_leaf_surface_conductance = leaf.calc_stomatal_conductance(
            residual_stomatal_conductance=residual_stomatal_conductance,
            maximum_stomatal_conductance=maximum_stomatal_conductance,
            absorbed_irradiance=absorbed_irradiance,
            shape_parameter=shape_parameter,
            stomatal_sensibility_to_water_status=stomatal_sensibility_to_water_status)

        leaf_fraction = calc_leaf_fraction(leaves_category,
                                           cumulative_leaf_area_index,
                                           direct_black_extinction_coefficient)

        leaf_layer_surface_conductance += (lumped_leaf_surface_conductance * leaf_fraction * sublayer_thickness)

    return leaf_layer_surface_conductance


def calc_leaf_layer_surface_resistance_to_vapor(leaves_category: str,
                                                incident_direct_irradiance: float,
                                                incident_diffuse_irradiance: float,
                                                upper_cumulative_leaf_area_index: float,
                                                lower_cumulative_leaf_area_index: float,
                                                stomatal_sensibility_to_water_status: float,
                                                leaf_scattering_coefficient: float,
                                                canopy_reflectance_to_direct_irradiance: float,
                                                canopy_reflectance_to_diffuse_irradiance: float,
                                                direct_extinction_coefficient: float,
                                                direct_black_extinction_coefficient: float,
                                                diffuse_extinction_coefficient: float,
                                                maximum_stomatal_conductance: float,
                                                residual_stomatal_conductance: float,
                                                shape_parameter: float,
                                                sublayers_number: int,
                                                stomatal_density_factor: int) -> float:
    """Calculates the bulk surface resistance of a leaf layer.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        incident_direct_irradiance: [W m-2ground] incident direct photosynthetically active radiation above the canopy
        incident_diffuse_irradiance: [W m-2ground] incident diffuse irradiance at the top of the canopy
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index above the considered layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index below the considered layer
        stomatal_sensibility_to_water_status: [-] stomatal closure fraction due to water stress
        leaf_scattering_coefficient: [-] leaf scattering coefficient
        canopy_reflectance_to_direct_irradiance: [-] canopy reflectance to direct (beam) irradiance
        canopy_reflectance_to_diffuse_irradiance: [-] canopy reflectance to diffuse irradiance
        direct_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance
        direct_black_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam)
            irradiance for black leaves
        diffuse_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of diffuse irradiance
        maximum_stomatal_conductance: [m h-1] maximum stomatal conductance
        residual_stomatal_conductance: [m h-1] residual (minimum) stomatal conductance
        shape_parameter: [W m-2leaf] an empirical parameter to regulate the shape of stomatal conductance response
            to absorbed photosynthetically active radiation
        sublayers_number: number of sublayers that are used to perform the numerical integral of the leaf-layer surface
            conductance equation
        stomatal_density_factor (int): [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    Returns:
        (float): [h m-1] bulk surface resistance of the leaf layer
    """

    args = {k: v for k, v in locals().items() if k != 'stomatal_density_factor'}
    surface_conductance = max(1.e-6, calc_leaf_layer_surface_conductance_to_vapor(**args))

    return stomatal_density_factor * (1.0 / surface_conductance)


def calc_leaf_fraction_per_leaf_layer(leaves_category: str,
                                      upper_cumulative_leaf_area_index: float,
                                      lower_cumulative_leaf_area_index: float,
                                      direct_black_extinction_coefficient: float) -> float:
    """Calculates the fraction of sunlit or shaded leaves at a given depth inside the canopy.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        direct_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of direct (beam) irradiance
            through a canopy of black leaves

    Returns:
        [-] fraction of sunlit or shaded leaves of the considered layer
    """
    leaf_layer_thickness = lower_cumulative_leaf_area_index - upper_cumulative_leaf_area_index
    sunlit_fraction_per_leaf_layer = calc_sunlit_fraction_per_leaf_layer(
        upper_cumulative_leaf_area_index, leaf_layer_thickness, direct_black_extinction_coefficient)
    if leaves_category == 'sunlit':
        return sunlit_fraction_per_leaf_layer
    else:
        return 1 - sunlit_fraction_per_leaf_layer


def calc_leaf_layer_net_longwave_radiation(leaves_category: str,
                                           canopy_top_net_longwave_radiation: float,
                                           upper_cumulative_leaf_area_index: float,
                                           lower_cumulative_leaf_area_index: float,
                                           direct_black_extinction_coefficient: float,
                                           diffuse_black_extinction_coefficient: float) -> float:
    """Calculates net long wave radiation exchange of a sunlit or shaded leaf layer.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        canopy_top_net_longwave_radiation: [W m-2ground] net long wave radiation at the top of the canopy
        upper_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index above the considered layer
        lower_cumulative_leaf_area_index: [m2leaf m-2ground] cumulative leaf area index below the considered layer
        direct_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of direct (beam) irradiance
            through a canopy of black leaves
        diffuse_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of diffuse irradiance
            through a canopy of black leaves

    Returns:
        [W m-2ground]: net long wave radiation exchange of a sunlit or shaded leaf layer

    References:
        Leuning et al. 1995
            Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies.
            Plant, Cell and Environment 18, 1183 - 1200.
    """
    extinction_coefficient = direct_black_extinction_coefficient + diffuse_black_extinction_coefficient
    sunlit_scaling_factor = (diffuse_black_extinction_coefficient / extinction_coefficient) * (
            exp(-extinction_coefficient * upper_cumulative_leaf_area_index) -
            exp(-extinction_coefficient * lower_cumulative_leaf_area_index))
    sunlit_net_longwave_radiation = canopy_top_net_longwave_radiation * sunlit_scaling_factor
    if leaves_category == 'sunlit':
        return sunlit_net_longwave_radiation
    else:
        lumped_layer_net_longwave_radiation = lumped_leaves.calc_leaf_layer_net_longwave_radiation(
            canopy_top_net_longwave_radiation=canopy_top_net_longwave_radiation,
            upper_cumulative_leaf_area_index=upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=lower_cumulative_leaf_area_index,
            diffuse_black_extinction_coefficient=diffuse_black_extinction_coefficient)
        return lumped_layer_net_longwave_radiation - sunlit_net_longwave_radiation
