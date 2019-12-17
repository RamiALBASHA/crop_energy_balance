from json import load

from crop_energy_balance import utils, params
from crop_energy_balance.formalisms import common

constants = params.Constants()


class Inputs:
    def __init__(self, inputs_path):
        inputs = load(open(str(inputs_path), mode='r'), encoding='utf-8')

        self.measurement_height = inputs['measurement_height']
        """[m] height at which the meteorological variables are measured"""

        self.canopy_height = inputs['canopy_height']
        """[m] height of the canopy"""

        self.soil_saturation_ratio = inputs['soil_saturation_ratio']
        """[-] ratio of actual to potential volumetric water content in the soil"""

        self.air_temperature = utils.celsius_to_kelvin(inputs['air_temperature'], constants.absolute_zero)
        """[K] temperature of the air"""

        self.wind_speed = inputs['wind_speed']
        """[m h-1] wind speed"""

        self.wind_speed_at_canopy_height = common.calc_wind_speed_at_canopy_height(
            self.wind_speed, self.canopy_height, self.measurement_height)
        """[m h-1] wind speed at canopy height"""

        self._air_vapor_pressure = inputs['vapor_pressure']
        """[kPa] vapor pressure of the air"""

        self.vapor_pressure_deficit = inputs['vapor_pressure_deficit']
        """[kPa] vapor pressure deficit of the air"""

        self.incident_par = inputs['incident_par']
        """[W m-2ground] incident photosynthetically active radiation"""

        self.net_shortwave_radiation = inputs['absorbed_global_irradiance']
        """[W m-2ground] absorbed net shortwave (global) irradiance per leaf layer

        Notes:
            Dictionary keys are integers indicating the order of the leaf layers.
            The uppermost layer must have the highest number while the lowermost layer has the lowest number.
        """

        self.sky_emissivity = common.calc_atmospheric_emissivity(self._air_vapor_pressure, self.air_temperature)
        """[-] sky longwave radiation emissivity"""

        self.psychrometric_constant = common.calc_psychrometric_constant(
            inputs['atmospheric_pressure'], constants.air_specific_heat_capacity,
            constants.latent_heat_for_vaporization, constants.vapor_to_dry_air_molecular_weight)
        """[kPa K-1] psychrometric constant
        
        See Also:
            :func:`calc_psychrometric_constant`
        """


class LumpedInputs(Inputs):
    def __init__(self, inputs_path):
        Inputs.__init__(self, inputs_path)
        inputs = load(open(str(inputs_path), mode='r'), encoding='utf-8')
        inputs = self._fmt_inputs(inputs)

        self.leaf_layers = inputs['leaf_layers']
        """[m2leaf m-2ground] dictionary of leaf area index per layer

        Notes:
            Dictionary keys are integers indicating the order of the leaf layers.
            The uppermost layer must have the highest number while the lowermost layer has the lowest number.
        """
        self.components_keys = sorted(list(self.leaf_layers.keys()) + [-1])

    @staticmethod
    def _fmt_inputs(inputs: dict):
        inputs['leaf_layers'] = {int(key): value for key, value in
                                 inputs['layered_leaf_area_index'].items()}

        inputs['absorbed__irradiance'] = {int(key): value for key, value in
                                          inputs['absorbed_global_irradiance'].items()}

        return inputs


class SunLitShadedInputs(Inputs):
    def __init__(self, inputs_path):
        Inputs.__init__(self, inputs_path)
