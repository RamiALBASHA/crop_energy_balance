from math import exp

from crop_energy_balance.formalisms import lumped_leaves, leaf
from crop_energy_balance.formalisms.irradiance import calc_leaf_fraction, calc_absorbed_irradiance
from crop_energy_balance.utils import discretize_linearly


def calc_leaf_layer_boundary_conductance(leaves_category: str,
                                         wind_speed_at_canopy_height: float,
                                         upper_leaf_area_index: float,
                                         lower_leaf_area_index: float,
                                         direct_black_extinction_coefficient: float,
                                         wind_speed_extinction_coefficient: float = 0.5,
                                         characteristic_length: float = 0.01,
                                         shape_parameter: float = 0.01) -> float:
    """Calculates bulk layer boundary layer conductance (for both sides of leaves).

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        upper_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        direct_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of direct (beam) irradiance
            through a canopy of black leaves
        wind_speed_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter

    Returns:
        [m h-1] bulk layer boundary layer conductance (for both sides of leaves)
    """
    leaf_boundary_conductance = leaf.calc_leaf_boundary_conductance(wind_speed_at_canopy_height,
                                                                    characteristic_length,
                                                                    shape_parameter)
    lumped_extinction_coefficient = 0.5 * wind_speed_extinction_coefficient + direct_black_extinction_coefficient
    sunlit_layer_scaling_factor = 1.0 / lumped_extinction_coefficient * (
            exp(-lumped_extinction_coefficient * upper_leaf_area_index) -
            exp(-lumped_extinction_coefficient * lower_leaf_area_index))

    sunlit_layer_boundary_conductance = leaf_boundary_conductance * sunlit_layer_scaling_factor

    if leaves_category == 'sunlit':
        return sunlit_layer_boundary_conductance
    elif leaves_category == 'sunlit-shaded':
        lumped_layer_boundary_conductance = lumped_leaves.calc_leaf_layer_boundary_conductance(
            wind_speed_at_canopy_height,
            upper_leaf_area_index,
            lower_leaf_area_index,
            wind_speed_extinction_coefficient,
            characteristic_length,
            shape_parameter)
        return lumped_layer_boundary_conductance - sunlit_layer_boundary_conductance


def calc_leaf_layer_boundary_resistance_to_vapor(leaves_category: str,
                                                 wind_speed_at_canopy_height: float,
                                                 upper_leaf_area_index: float,
                                                 lower_leaf_area_index: float,
                                                 direct_black_extinction_coefficient: float,
                                                 wind_speed_extinction_coefficient: float,
                                                 characteristic_length: float,
                                                 shape_parameter: float,
                                                 stomatal_density_factor: int) -> float:
    """Calculates the bulk leaf layer resistance to water vapor transfer.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        upper_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the top of the layer
        lower_leaf_area_index: [m2leaf m-2ground] cumulative leaf layer index at the bottom of the layer
        direct_black_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of direct (beam) irradiance
            through a canopy of black leaves
        wind_speed_extinction_coefficient: [m2ground m-2leaf] extinction coefficient of wind speed inside the canopy
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter
        stomatal_density_factor (int): [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2

    Returns:
        [h m-1] bulk leaf layer resistance to water vapor transfer
    """
    args = {k: v for k, v in locals().items() if k != 'stomatal_density_factor'}
    boundary_conductance = calc_leaf_layer_boundary_conductance(**args)

    return 1.0 / (boundary_conductance / stomatal_density_factor)


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
    """Calculates the bulk surface conductance of a leaf layer.

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
        sublayers_number: number of sublayers that are used to performe numerical integral of the leaf-layer surface
            conductance equation

    Returns:
        [m h-1] bulk surface conductance of the leaf layer
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

        stomatal_sensibility_to_absorbed_irradiance = (absorbed_irradiance / (shape_parameter + absorbed_irradiance))

        lumped_leaf_surface_conductance = residual_stomatal_conductance + maximum_stomatal_conductance * (
                stomatal_sensibility_to_water_status * stomatal_sensibility_to_absorbed_irradiance)

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
                                                shape_parameter: float = 105,
                                                sublayers_number: int = 5) -> float:
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
        sublayers_number: number of sublayers that are used to performe numerical integral of the leaf-layer surface
            conductance equation

    Returns:
        (float): [h m-1] bulk surface resistance of the leaf layer
    """
    return 1.0 / max(1.e-6, calc_leaf_layer_surface_conductance_to_vapor(**locals()))


def calc_leaf_layer_boundary_resistance_to_heat(leaf_layer_boundary_conductance: float) -> float:
    """Calculates the bulk leaf layer resistance to heat transfer.

    Args:
        leaf_layer_boundary_conductance: [m h-1] Calculates bulk layer boundary layer conductance (for both blade sides)

    Returns:
        [h m-1] bulk leaf layer resistance to heat transfer

    """
    return 1.0 / max(1.e-6, leaf_layer_boundary_conductance)
