from json import load


class Params:
    def __init__(self, params_path):
        self.user_params = load(open(str(params_path), mode='r'), encoding='utf-8')

        self.simulation = Simulation(self.user_params)
        self.numerical_resolution = NumericalResolution(self.user_params)


class Constants:
    def __init__(self):
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


class Simulation:
    def __init__(self, data):
        self.vapor_pressure_deficit_coefficient = data['vpd_coeff']
        """[kPa] empirical parameter regulating the shape of stomata response to air vapor pressure deficit"""

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
        """ [-] 1 for amphistomatal leaves (stomata on both sides of the blade), otherwise 2"""

        self.extinction_coefficient = None
        """[m2ground m-2leaf] extinction coefficient of lumped direct and diffuse irradiance"""

        self.diffuse_black_extinction_coefficient = None
        """[m2ground m-2leaf] extinction coefficient of diffuse photosynthetically active radiation for black leaves"""


class LumpedSimulation(Simulation):
    def __init__(self, data):
        Simulation.__init__(self, data)

        self.global_extinction_coefficient = data['global_extinction_coefficient']
        """[m2ground m-2leaf] extinction coefficient of global irradiance"""

        self.diffuse_extinction_coef = data['diffuse_extinction_coef']
        """[m2ground m-2leaf] extinction coefficient of diffuse photosynthetically active radiation for black leaves"""


class SunlitShadedSimulation(Simulation):
    def __init__(self,
                 data: dict):
        Simulation.__init__(self, data)

        self.direct_black_extinction_coefficient = None

    def update(self, data):
        self.direct_black_extinction_coefficient = data['global_extinction_coefficient']
        """[m2ground m-2leaf] extinction coefficient of global irradiance"""


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
