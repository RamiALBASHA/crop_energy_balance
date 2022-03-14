from json import load
from pathlib import Path

from crop_irradiance.uniform_crops.formalisms.sunlit_shaded_leaves import (calc_canopy_reflectance_to_direct_irradiance,
                                                                           calc_diffuse_extinction_coefficient,
                                                                           calc_direct_extinction_coefficient,
                                                                           calc_direct_black_extinction_coefficient)

from crop_energy_balance.formalisms.weather import calc_atmospheric_emissivity


class Params:
    def __init__(self,
                 params_dict: dict = None,
                 params_path: Path = None):
        if params_dict:
            self._user_params = params_dict
        else:
            with open(str(params_path), mode='r', encoding='utf-8') as f:
                self._user_params = load(f)

        self.simulation = Simulation(self._user_params)

        self.numerical_resolution = NumericalResolution(self._user_params)

    def update(self,
               inputs):
        self.simulation.update(inputs=inputs)


class Constants:
    def __init__(self):
        self.gravitational_acceleration = 9.81 * (3600 ** 2)
        """[m h-2] gravitational acceleration.
        """

        self.von_karman = 0.41
        """[-] von Karman constant"""

        self.stefan_boltzmann = 5.67e-8
        """[W m-2 K-4] Stefan-Boltzmann constant"""

        self.absolute_zero = -273.15
        """[°C] temeprature at absolute zero"""

        self.latent_heat_for_vaporization = 0.678
        """[W h g-1] latent heat for vaporization"""

        self.psychrometric_constant = 0.066
        """[kPa K-1] psychrometric constant"""

        self.air_specific_heat_capacity = 2.8e-4
        """[W h g-1 K-1] specific heat capacity of the air under a constant pressure

        References:
            Allen et al. 1998
                FAO Irrigation and Drainage Paper No. 56.
                Eq. 8
        """

        self.vapor_to_dry_air_molecular_weight = 0.622
        """[-] ratio of the molecular weights of water vapor to dry air

        References:
            Allen et al. 1998
                FAO Irrigation and Drainage Paper No. 56.
                Eq. 8
        """

        self.air_density = 1185.0
        """[g m-3] dry air density
        """

        self.ideal_gas_constant = 8.2057 * 1.e-5
        """[m3 atm mol−1 K−1] Ideal gas constant
        """

        self.molecular_diffusivity_water_vapor = 3600. * 1e-6 * 24.9
        """[m2 h-1] Molucular diffusivity for water vapor at 25 degrees C.

        References:
            Monteith and Unsworth (2013).
                Principles of Environmental Physics (Fourth Edition)
                Academic Press, pp 289 - 320
                Table A.3
        """


class Simulation:
    def __init__(self, data):
        self.stomatal_sensibility = data['stomatal_sensibility']
        """A dictionary of stomatal sensibility parameters (multiple models are handled):
            - key: name of the model (e.g. 'leuning', 'tuzet', 'misson')
            - value: dictionary :
                - key: parameter name
                - value: parameter value
        """

        self.soil_aerodynamic_resistance_shape_parameter = data['soil_aerodynamic_resistance_shape_parameter']
        """[-] soil aerodynamic resistance shape parameter.
        References:
            Choudhury, Monteith, 1988
                A four-layer model for the heat budget of homogeneous land surfaces.
                Quarterly Journal of the Royal Meteorological Society 114, 373 – 398
        """

        self.soil_roughness_length_for_momentum = data['soil_roughness_length_for_momentum']
        """[m] soil roughness length for momentum.

        References:
            Choudhury, Monteith, 1988
                A four-layer model for the heat budget of homogeneous land surfaces.
                Quarterly Journal of the Royal Meteorological Society 114, 373 – 398
        """

        self.leaf_characteristic_length = data['leaf_characteristic_length']
        """[m] leaf characteristic length in the direction of wind
        
        See Also:
            Yin and van Laar, 1994
                Crop Systems Dynamics: An ecophysiological simulation model for genotype-by-environment interactions.
                Wageningen Academic Publishers pp. 169
        """

        self.leaf_boundary_layer_shape_parameter = data['leaf_boundary_layer_shape_parameter']
        """[m s-1/2] empirical shape parameter for the calculation of leaf-level boundary-layer conductance.

        See Also:
            Yin and van Laar, 1994
                Crop Systems Dynamics: An ecophysiological simulation model for genotype-by-environment interactions.
                Wageningen Academic Publishers pp. 169
        """

        self.wind_speed_extinction_coef = data['wind_speed_extinction_coef']
        """[m2ground m-2leaf] extinction coefficient of wind speed through the canopy.

        See Also:
            Inoue, 1963
                On the turbulent structure of air flow within crop canopies.
                Journal of Meteorological Society of Japan 41, 317 – 326
        """

        self.maximum_stomatal_conductance = data['maximum_stomatal_conductance']
        """[m h-1] maximum stomatal conductance.

        Notes:
            resistance [m2 s-1 mol-1] = 0.024 * 3600 / resistance [h m-1],
                cf. Eq. 3.15 in Monteith and Unsworth (2013). Principles of Environmental Physics, 4th Edition. 
        """

        self.residual_stomatal_conductance = data['residual_stomatal_conductance']
        """[m h-1] maximum stomatal conductance.
        Notes:
            resistance [m2 s-1 mol-1] = 0.024 * 3600 / resistance [h m-1],
                cf. Eq. 3.15 in Monteith and Unsworth (2013). Principles of Environmental Physics, 4th Edition. 
        """

        self.leaf_emissivity = data['leaf_emissivity']
        """[-] leaf emissivity to longwave energy"""

        self.soil_emissivity = data['soil_emissivity']
        """[-] soil emissivity to longwave energy"""

        self.absorbed_par_50 = data['absorbed_par_50']
        """[W m-2leaf] absorbed photosynthetically active radiation flux density at which stomatal conductance is half
        its maximum value.

        See Also:
            Uddling and Pleijel, 2006.
                Changes in stomatal conductance and net photosynthesis during phenological development in spring wheat:
                    implications for gas exchange modelling
                International Journal of Biometeorology 51, 37 – 48
                Extracted from Figure 2a.
        """

        self.soil_resistance_to_vapor_shape_parameter_1 = data['soil_resistance_to_vapor_shape_parameter_1']
        """[-] empirical parameter to control the shape of the relationship between soil resistance to water vapor
        transfer and soil saturation fraction.
        """

        self.soil_resistance_to_vapor_shape_parameter_2 = data['soil_resistance_to_vapor_shape_parameter_2']
        """[-] empirical parameter to control the shape of the relationship between soil resistance to water vapor
        transfer and soil saturation fraction.
        """

        self.stomatal_density_factor = data['stomatal_density_factor']
        """[-] 1 for amphistomatal leaves (stomata on both sides of the blade) or 2 for hypostomatal leaves
        (stomata on one side of the blade)."""

        self.leaf_scattering_coefficient = data['leaf_scattering_coefficient']
        """[-] leaf scattering coefficient"""

        self.canopy_reflectance_to_diffuse_irradiance = 0.057
        """[-] canopy reflectance to diffuse irradiance"""

        if 'leaf_angle_distribution_factor' in data.keys():
            self.leaf_angle_distribution_factor = data['leaf_angle_distribution_factor']
        else:
            self.leaf_angle_distribution_factor = 0.9773843811168246
        """[-] factor describing leaf angle distribution (for spherical distributions its value equals
        rad(56) = 0.9773843811168246"""

        if 'clumping_factor' in data.keys():
            self.clumping_factor = data['clumping_factor']
        else:
            self.clumping_factor = 1
        """[-] clumping factor to describe the spatial dependency of the positions of the leaves"""

        self.sublayers_number = 100
        """[-] number of sublayers that are used to perform the numerical integral of the leaf-layer surface conductance
        equation
        """

        self.canopy_reflectance_to_direct_irradiance = None
        """[-] canopy reflectance to direct (beam) irradiance"""

        self.direct_extinction_coefficient = None
        """[m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance"""

        self.direct_black_extinction_coefficient = None
        """[m2ground m-2leaf] the extinction coefficient of direct (beam) irradiance for black leaves"""

        self.diffuse_extinction_coefficient = None
        """[m2ground m-2leaf] the extinction coefficient of diffuse irradiance"""

        self.diffuse_black_extinction_coefficient = None
        """[m2ground m-2leaf] extinction coefficient of diffuse photosynthetically active radiation for black leaves"""

        self.drag_coefficient = 0.2
        """[m2ground m-2leaf] drag coefficient"""

        self.ratio_heat_to_momentum_canopy_roughness_lengths = 1 / 7.4
        """[-] Ratio of canopy's heat to momentum roughness lengths.
        Indicative values are:
            * 1/10 for reference grass crop (Shuttleworth, 2007. Hydrol. Earth Syst. Sci. 11, 210 - 244)
            * 1/7.4 for wheat (Kimball et al., 2015. Climatology and Water Management 107, 129 - 141)
        """

        self.richardon_threshold_free_convection = -0.8
        """[-] Richardson number threshold below which flux is assumed to occur under free convection.
        Note:
        Indicative values are:
            * -1.0 (Monteith and Unsworth, 2004. cf. description belw Eq. 16.45)
            * -0.8 (CanopyT, Webber et al., 2016)

        """

        self.atmospheric_emissivity = None
        """[-] sky longwave radiation emissivity"""

        if 'atmospheric_emissivity_model' in data.keys():
            self.atmospheric_emissivity_model = data['atmospheric_emissivity_model']
        else:
            self.atmospheric_emissivity_model = 'brutsaert_1975'
        """Name of the model to be used for calculating sky longwave radiation emissivity"""

    def update(self,
               inputs):
        self.direct_extinction_coefficient = calc_direct_extinction_coefficient(
            solar_inclination=inputs.solar_inclination,
            leaf_scattering_coefficient=self.leaf_scattering_coefficient,
            leaf_angle_distribution_factor=self.leaf_angle_distribution_factor,
            clumping_factor=self.clumping_factor)

        self.direct_black_extinction_coefficient = calc_direct_black_extinction_coefficient(
            solar_inclination=inputs.solar_inclination,
            leaf_angle_distribution_factor=self.leaf_angle_distribution_factor,
            clumping_factor=self.clumping_factor)

        self.diffuse_extinction_coefficient, self.diffuse_black_extinction_coefficient = (
            calc_diffuse_extinction_coefficient(
                leaf_area_index=sum(inputs.leaf_layers.values()),
                leaf_angle_distribution_factor=self.leaf_angle_distribution_factor,
                clumping_factor=self.clumping_factor,
                leaf_scattering_coefficient=self.leaf_scattering_coefficient,
                sky_sectors_number=3,
                sky_type='soc'))

        self.canopy_reflectance_to_direct_irradiance = calc_canopy_reflectance_to_direct_irradiance(
            direct_black_extinction_coefficient=self.direct_black_extinction_coefficient,
            leaf_scattering_coefficient=self.leaf_scattering_coefficient)

        self.atmospheric_emissivity = calc_atmospheric_emissivity(
            model=self.atmospheric_emissivity_model,
            air_vapor_pressure=inputs.air_vapor_pressure,
            air_temperature=inputs.air_temperature)
        """[-] sky longwave radiation emissivity"""


class NumericalResolution:
    def __init__(self, data):
        self.step_fraction = data['step_fraction']
        """[-] fraction of the entire temperature step (`actual_value - previous_value`) to be used"""

        self.acceptable_temperature_error = data['acceptable_temperature_error']
        """[°C] acceptable temperature estimation error between two consecutive iterations.
        
        See Also:
            Maes and Steppe, 2012
                Estimating evapotranspiration and drought stress with ground-based thermal remote sensing in
                    agriculture: a review.
                Journal of Experimental Botany 63, 4671 – 4712
        """

        self.maximum_iteration_number = data['maximum_iteration_number']
        """[-] maximum number of iterations to solve the energy budget"""
