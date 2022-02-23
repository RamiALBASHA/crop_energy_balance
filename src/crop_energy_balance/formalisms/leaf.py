from math import exp


def calc_forced_convection_condutance(wind_speed_at_canopy_height: float,
                                      characteristic_length: float = 0.01,
                                      shape_parameter: float = 0.01) -> float:
    """Calculates boundary layer conductance under forced convection (for both sides of leaves) at the scale of an
    individual leaf.

    Args:
        wind_speed_at_canopy_height: [m s-1] local wind speed in the vicinity of the leaf
        characteristic_length: [m] characteristic leaf length in the direction of the wind
        shape_parameter: [m s-0.5] an empirical shape parameter

    Returns:
        [m h-1] bulk boundary layer conductance uner forced convection (for both sides of leaves) at the scale of an
        individual leaf

    """
    return 3600 * shape_parameter * (wind_speed_at_canopy_height / characteristic_length) ** 0.5


def calc_stomatal_sensibility(model_args: dict) -> float:
    """Calculates the effect of water status on stomatal aperture following one or several models.

    Args:
        model_args: A dictionary of stomatal sensibility parameters (multiple models are handled)
            - key: name of the model (one of 'leuning', 'tuzet', 'misson')
            - value: dictionary :
                - key: parameter name
                - value: parameter value

    Returns:
        [-] reduction fraction of the maximum stomatal aperture

    Note:
        If multiple models are provided, the resulting reduction fraction is the product of all reduction fractions
            from all models

    """
    reduction_factor = 1
    for name, kwargs in model_args.items():
        if name == 'leuning':
            reduction_factor *= calc_stomatal_sensibility_leuning(**kwargs)
        elif name == 'tuzet':
            reduction_factor *= calc_stomatal_sensibility_tuzet(**kwargs)
        elif name == 'misson':
            reduction_factor *= calc_stomatal_sensibility_misson(**kwargs)
        else:
            raise ValueError(f'Model "name" is not recognized.')

    return reduction_factor


def calc_stomatal_sensibility_leuning(vapor_pressure_deficit: float,
                                      d_0: float) -> float:
    """Calculates the effect of air vapor pressure deficit on stomatal aperture following Leuning et al. (1995).

    Args:
        vapor_pressure_deficit: [kPa] vapor pressure deficit of the air at canopy source height
        d_0: [kPa] empirical shape parameter

    Returns:
        [-] reduction fraction of the maximum stomatal aperture

    Reference:
        Leuning (1995)
            A critical appraisal of a combined stomatal photosynthesis model for C3 plants.
            Plant, Cell and Environment 18, 339355.

    """
    assert (d_0 != 0), 'The value of `shape_parameter` must be greater than zero.'
    return 1.0 / (1.0 + vapor_pressure_deficit / d_0)


def calc_stomatal_sensibility_tuzet(water_potential: float, psi_ref: float, steepness: float) -> float:
    """Calculates the effect of leaf water potential on stomatal aperture following Tuzet et al. (2003).

    Args:
        water_potential: [MPa] bulk leaf water potential
        psi_ref: [MPa] empirical shape parameter
        steepness: [MPa-1] steepness of the sigmoidal function

    Returns:
        [-] reduction fraction of the maximum stomatal aperture

    Reference:
        Tuzet et al. (2003)
            A coupled model of stomatal conductance, photosynthesis and transpiration.
            Plant, Cell and Environment 26, 10971116.

    """
    return (1. + exp(steepness * psi_ref)) / (1. + exp(float(steepness * (psi_ref - water_potential))))


def calc_stomatal_sensibility_misson(water_potential: float, psi_half_aperture: float, steepness: float) -> float:
    """Calculates the effect of soil water potential on stomatal aperture following Misson et al. (2004).

    Args:
        water_potential: [MPa] bulk leaf water potential
        psi_half_aperture: [MPa] value of psi at which the maximum stomatal conductance is reduced by 50 %
        steepness: [-] steepness of the sigmoidal function

    Returns:
        [-] reduction fraction of the maximum stomatal aperture

    Reference:
        Misson et al. (2004)
            A comparison of three approaches to modeling leaf gas exchange in annually drought-stressed ponderosa pine
                forest.
            Tree Physiology 24, 529 541.

    """
    return 1. / (1. + (water_potential / psi_half_aperture) ** steepness)


def calc_stomatal_conductance(residual_stomatal_conductance: float,
                              maximum_stomatal_conductance: float,
                              absorbed_irradiance: float,
                              shape_parameter: float,
                              stomatal_sensibility_to_water_status: float) -> float:
    """Calculates stomatal conductivity to water vapor for a unit leaf surface area.

    Args:
        residual_stomatal_conductance: [m h-1] residual (minimum) stomatal conductance
        maximum_stomatal_conductance: [m h-1] maximum stomatal conductance
        absorbed_irradiance: [W m-2leaf] absorbed photosynthetically active radiation (PAR)
        shape_parameter: [W m-2leaf] an empirical parameter that regulates the response of stomatal conductance to the
            absorbed PAR
        stomatal_sensibility_to_water_status: [-] stomatal closure fraction due to water stress

    Returns:
        [m h-1] stomatal conductivity to water vapor for a unit leaf surface area
    """
    stomatal_sensibility_to_absorbed_irradiance = (absorbed_irradiance / (shape_parameter + absorbed_irradiance))
    return residual_stomatal_conductance + maximum_stomatal_conductance * (
            stomatal_sensibility_to_absorbed_irradiance * stomatal_sensibility_to_water_status)
