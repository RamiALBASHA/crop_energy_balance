from crop_irradiance.uniform_crops.formalisms.sunlit_shaded_leaves import (
    calc_sunlit_fraction,
    calc_shaded_fraction,
    calc_absorbed_direct_irradiance,
    calc_absorbed_diffuse_irradiance_at_given_depth,
    calc_absorbed_scattered_irradiance_at_given_depth)


def calc_absorbed_irradiance_by_sunlit_leaves_per_leaf_layer_at_given_depth(incident_direct_irradiance: float,
                                                                            incident_diffuse_irradiance: float,
                                                                            cumulative_leaf_area_index: float,
                                                                            leaf_scattering_coefficient: float,
                                                                            canopy_reflectance_to_direct_irradiance: float,
                                                                            canopy_reflectance_to_diffuse_irradiance: float,
                                                                            direct_extinction_coefficient: float,
                                                                            direct_black_extinction_coefficient: float,
                                                                            diffuse_extinction_coefficient: float) -> float:
    """Calculates the absorbed irradiance by sunlit leaves at a given depth, per unit leaf area (depth-independent).

    Args:
        incident_direct_irradiance: [W m-2ground] incident direct (beam) irradiance at the top of the canopy
        incident_diffuse_irradiance: [W m-2ground] incident diffuse irradiance at the top of the canopy
        cumulative_leaf_area_index: [m2leaf m-2ground] cumulative downwards leaf area index at the top of the
            considered layer
        leaf_scattering_coefficient: [-] leaf scattering coefficient
        canopy_reflectance_to_direct_irradiance: [-] canopy reflectance to direct (beam) irradiance
        canopy_reflectance_to_diffuse_irradiance: [-] canopy reflectance to diffuse irradiance
        direct_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance
        direct_black_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam)
            irradiance for black leaves
        diffuse_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of diffuse irradiance

    Returns:
        [W m-2leaf] the absorbed irradiance by sunlit leaves at a given depth, per unit leaf area
    """
    absorbed_direct_irradiance = calc_absorbed_direct_irradiance(
        incident_direct_irradiance,
        leaf_scattering_coefficient,
        direct_black_extinction_coefficient)
    absorbed_diffuse_irradiance = calc_absorbed_diffuse_irradiance_at_given_depth(
        incident_diffuse_irradiance,
        cumulative_leaf_area_index,
        canopy_reflectance_to_diffuse_irradiance,
        diffuse_extinction_coefficient)
    absorbed_scattered_irradiance = calc_absorbed_scattered_irradiance_at_given_depth(
        incident_direct_irradiance,
        cumulative_leaf_area_index,
        direct_extinction_coefficient,
        direct_black_extinction_coefficient,
        canopy_reflectance_to_direct_irradiance,
        leaf_scattering_coefficient)

    return absorbed_direct_irradiance + absorbed_diffuse_irradiance + absorbed_scattered_irradiance


def calc_absorbed_irradiance_by_shaded_leaves_per_leaf_layer_at_given_depth(incident_direct_irradiance: float,
                                                                            incident_diffuse_irradiance: float,
                                                                            cumulative_leaf_area_index: float,
                                                                            leaf_scattering_coefficient: float,
                                                                            canopy_reflectance_to_direct_irradiance: float,
                                                                            canopy_reflectance_to_diffuse_irradiance: float,
                                                                            direct_extinction_coefficient: float,
                                                                            direct_black_extinction_coefficient: float,
                                                                            diffuse_extinction_coefficient: float) -> float:
    """Calculates the absorbed irradiance by shaded leaves at a given depth, per unit leaf area (depth-independent).

    Args:
        incident_direct_irradiance: [W m-2ground] incident direct (beam) irradiance at the top of the canopy
        incident_diffuse_irradiance: [W m-2ground] incident diffuse irradiance at the top of the canopy
        cumulative_leaf_area_index: [m2leaf m-2ground] cumulative downwards leaf area index at the top of the
            considered layer
        leaf_scattering_coefficient: [-] leaf scattering coefficient
        canopy_reflectance_to_direct_irradiance: [-] canopy reflectance to direct (beam) irradiance
        canopy_reflectance_to_diffuse_irradiance: [-] canopy reflectance to diffuse irradiance
        direct_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance
        direct_black_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam)
            irradiance for black leaves
        diffuse_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of diffuse irradiance

    Returns:
        [W m-2leaf] the absorbed irradiance by shaded leaves at a given depth, per unit leaf area
    """
    absorbed_diffuse_irradiance = calc_absorbed_diffuse_irradiance_at_given_depth(
        incident_diffuse_irradiance,
        cumulative_leaf_area_index,
        canopy_reflectance_to_diffuse_irradiance,
        diffuse_extinction_coefficient)
    absorbed_scattered_irradiance = calc_absorbed_scattered_irradiance_at_given_depth(
        incident_direct_irradiance,
        cumulative_leaf_area_index,
        direct_extinction_coefficient,
        direct_black_extinction_coefficient,
        canopy_reflectance_to_direct_irradiance,
        leaf_scattering_coefficient)

    return absorbed_diffuse_irradiance + absorbed_scattered_irradiance


def calc_absorbed_irradiance(leaves_category: str,
                             incident_direct_irradiance: float,
                             incident_diffuse_irradiance: float,
                             cumulative_leaf_area_index: float,
                             leaf_scattering_coefficient: float,
                             canopy_reflectance_to_direct_irradiance: float,
                             canopy_reflectance_to_diffuse_irradiance: float,
                             direct_extinction_coefficient: float,
                             direct_black_extinction_coefficient: float,
                             diffuse_extinction_coefficient: float) -> float:
    """Calculates the absorbed irradiance by sunlit or shaded leaves at a given depth, per unit leaf area
        (depth-independent).

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        incident_direct_irradiance: [W m-2ground] incident direct (beam) irradiance at the top of the canopy
        incident_diffuse_irradiance: [W m-2ground] incident diffuse irradiance at the top of the canopy
        cumulative_leaf_area_index: [m2leaf m-2ground] cumulative downwards leaf area index at the top of the
            considered layer
        leaf_scattering_coefficient: [-] leaf scattering coefficient
        canopy_reflectance_to_direct_irradiance: [-] canopy reflectance to direct (beam) irradiance
        canopy_reflectance_to_diffuse_irradiance: [-] canopy reflectance to diffuse irradiance
        direct_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance
        direct_black_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam)
            irradiance for black leaves
        diffuse_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of diffuse irradiance

    Returns:
        [W m-2leaf] the absorbed irradiance by sunlit or shaded leaves at a given depth, per unit leaf area
    """
    if leaves_category == 'sunlit':
        return calc_absorbed_irradiance_by_sunlit_leaves_per_leaf_layer_at_given_depth(
            incident_direct_irradiance,
            incident_diffuse_irradiance,
            cumulative_leaf_area_index,
            leaf_scattering_coefficient,
            canopy_reflectance_to_direct_irradiance,
            canopy_reflectance_to_diffuse_irradiance,
            direct_extinction_coefficient,
            direct_black_extinction_coefficient,
            diffuse_extinction_coefficient)
    elif leaves_category == 'shaded':
        return calc_absorbed_irradiance_by_shaded_leaves_per_leaf_layer_at_given_depth(
            incident_direct_irradiance,
            incident_diffuse_irradiance,
            cumulative_leaf_area_index,
            leaf_scattering_coefficient,
            canopy_reflectance_to_direct_irradiance,
            canopy_reflectance_to_diffuse_irradiance,
            direct_extinction_coefficient,
            direct_black_extinction_coefficient,
            diffuse_extinction_coefficient)


def calc_leaf_fraction(leaves_category: str,
                       cumulative_leaf_area_index: float,
                       direct_black_extinction_coefficient: float) -> float:
    """Calculates the fraction of sunlit or shaded leaves at a given depth inside the canopy.

    Args:
        leaves_category: one of ('sunlit', 'shaded')
        cumulative_leaf_area_index: [m2leaf m-2ground] cumulative downwards leaf area index
        direct_black_extinction_coefficient: [m2ground m-2leaf] the extinction coefficient of direct (beam)
            irradiance for black leaves

    Returns:
        [-] fraction of sunlit or shaded leaves at a given depth inside the canopy
    """

    if leaves_category == 'sunlit':
        return calc_sunlit_fraction(cumulative_leaf_area_index, direct_black_extinction_coefficient)
    elif leaves_category == 'shaded':
        return calc_shaded_fraction(cumulative_leaf_area_index, direct_black_extinction_coefficient)
