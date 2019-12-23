from json import load

from crop_energy_balance.formalisms import canopy, weather

from crop_energy_balance import utils, params

constants = params.Constants()


class Inputs:
    def __init__(self, inputs_path):
        with open(str(inputs_path), mode='r') as f:
            inputs = load(f, encoding='utf-8')

        inputs = self._fmt_inputs(inputs)

        self.measurement_height = inputs['measurement_height']
        """[m] height at which the meteorological variables are measured"""

        self.canopy_height = inputs['canopy_height']
        """[m] height of the canopy"""

        self.soil_saturation_ratio = inputs['soil_saturation_ratio']
        """[-] ratio of actual to potential volumetric water content in the soil"""

        self.air_temperature = utils.convert_celsius_to_kelvin(inputs['air_temperature'], constants.absolute_zero)
        """[K] temperature of the air"""

        self.wind_speed = inputs['wind_speed']
        """[m h-1] wind speed"""

        self.wind_speed_at_canopy_height = canopy.calc_wind_speed_at_canopy_height(
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

        self.atmospheric_emissivity = weather.calc_atmospheric_emissivity(self._air_vapor_pressure,
                                                                          self.air_temperature)
        """[-] sky longwave radiation emissivity"""

        self.psychrometric_constant = weather.calc_psychrometric_constant(
            inputs['atmospheric_pressure'], constants.air_specific_heat_capacity,
            constants.latent_heat_for_vaporization, constants.vapor_to_dry_air_molecular_weight)
        """[kPa K-1] psychrometric constant
        
        See Also:
            :func:`calc_psychrometric_constant`
        """
        self.leaf_layers = inputs['leaf_layers']
        """[m2leaf m-2ground] dictionary of leaf area index per layer

        Notes:
            Dictionary keys are integers indicating the order of the leaf layers.
            The uppermost layer must have the highest number while the lowermost layer has the lowest number.
        """
        self.absorbed_irradiance = self._calc_absorbed_irradiance(inputs)
        """[W_{PAR} m-2ground] dictionary of absobed photosynthetically active radiation per layer layer

        Notes:
            For lumped leaves, the absorbed irradiance per leaf layer is a dictionary having the key 'lumped'
            For sunlit and shaded leaves, the absorbed irradiance per leaf layer is a dictionary having the keys
                'sunlit' and 'shaded'
            Dictionary keys are integers indicating the order of the leaf layers.
            The uppermost layer must have the highest number while the lowermost layer has the lowest number.
        """

        self.components_keys = sorted(list(self.leaf_layers.keys()) + [-1])

    @staticmethod
    def _calc_absorbed_irradiance(inputs) -> dict:
        res = {}
        for component_key, global_irradiance in inputs['absorbed_global_irradiance'].items():
            res[component_key] = {
                leaves_category: utils.convert_global_irradiance_into_photosynthetically_active_radition(value)
                for leaves_category, value in global_irradiance.items()}
        return res

    @staticmethod
    def _fmt_inputs(inputs: dict):
        inputs['leaf_layers'] = {int(key): value for key, value in
                                 inputs['leaf_layers'].items()}
        inputs['absorbed_global_irradiance'] = {int(key): value for key, value in
                                                inputs['absorbed_global_irradiance'].items()}
        return inputs
