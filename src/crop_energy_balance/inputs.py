from pathlib import Path
from json import load

from crop_energy_balance.formalisms import canopy, weather

from crop_energy_balance import utils, params

constants = params.Constants()


class Inputs:
    def __init__(self, inputs_path: Path):
        with open(str(inputs_path), mode='r') as f:
            inputs = load(f, encoding='utf-8')

        self._inputs = self._fmt_inputs(inputs)

        self.measurement_height = self._inputs['measurement_height']
        """[m] height at which the meteorological variables are measured"""

        self.canopy_height = self._inputs['canopy_height']
        """[m] height of the canopy"""

        self.soil_saturation_ratio = self._inputs['soil_saturation_ratio']
        """[-] ratio of actual to potential volumetric water content in the soil"""

        self.air_temperature = utils.convert_celsius_to_kelvin(self._inputs['air_temperature'], constants.absolute_zero)
        """[K] temperature of the air"""

        self.wind_speed = self._inputs['wind_speed']
        """[m h-1] wind speed"""

        self.wind_speed_at_canopy_height = canopy.calc_wind_speed_at_canopy_height(
            self.wind_speed, self.canopy_height, self.measurement_height)
        """[m h-1] wind speed at canopy height"""

        self._air_vapor_pressure = self._inputs['vapor_pressure']
        """[kPa] vapor pressure of the air"""

        self.vapor_pressure_deficit = self._inputs['vapor_pressure_deficit']
        """[kPa] vapor pressure deficit of the air"""

        self.atmospheric_emissivity = weather.calc_atmospheric_emissivity(self._air_vapor_pressure,
                                                                          self.air_temperature)
        """[-] sky longwave radiation emissivity"""

        self.psychrometric_constant = weather.calc_psychrometric_constant(
            self._inputs['atmospheric_pressure'], constants.air_specific_heat_capacity,
            constants.latent_heat_for_vaporization, constants.vapor_to_dry_air_molecular_weight)
        """[kPa K-1] psychrometric constant

        See Also:
            :func:`calc_psychrometric_constant`
        """

        self.leaf_layers = self._inputs['leaf_layers']
        """[m2leaf m-2ground] dictionary of leaf area index per layer

        Notes:
            Dictionary keys are integers indicating the order of the leaf layers.
            The uppermost layer must have the highest number while the lowermost layer has the lowest number.
        """

        self.absorbed_irradiance = self._inputs['absorbed_photosynthetically_active_radition']
        """[W_{PAR} m-2ground] dictionary of absorbed photosynthetically active radiation per crop component

        Notes:
            For lumped leaves, the absorbed irradiance per leaf layer is a dictionary having the key 'lumped'
            For sunlit and shaded leaves, the absorbed irradiance per leaf layer is a dictionary having the keys
                'sunlit' and 'shaded'
            Dictionary keys are integers indicating the order of the leaf layers.
            The uppermost layer must have the highest number while the lowermost layer has the lowest number.
        """

        self.net_shortwave_radiation = self._calc_net_shortwave_irradiance()
        """[W m-2ground] absorbed net shortwave (global) irradiance per crop component

        Notes:
            Dictionary keys are integers indicating the order of the crop components.
            The uppermost component must have the highest number while the lowermost component, i.e. soil, has the
                lowest number which must be equal to -1.
        """

        self.components_keys = sorted(list(self.leaf_layers.keys()) + [-1])

    def _calc_net_shortwave_irradiance(self) -> dict:
        res = {}
        for component_key, photosynthetically_active_radiation in self.absorbed_irradiance.items():
            for leaves_category, value in photosynthetically_active_radiation.items():
                res[component_key] = {
                    leaves_category: utils.convert_photosynthetically_active_radition_into_global_radiation(value)}
        return res

    @staticmethod
    def _fmt_inputs(inputs: dict):
        for k in ('leaf_layers', 'absorbed_photosynthetically_active_radition'):
            inputs[k] = {int(key): value for key, value in inputs[k].items()}
        return inputs


class LumpedInputs(Inputs):
    def __init__(self, inputs_path: Path):
        Inputs.__init__(self, inputs_path)


class SunlitShadedInputs(Inputs):
    def __init__(self, inputs_path: Path):
        Inputs.__init__(self, inputs_path)

        self.incident_par = self._inputs['incident_photosynthetically_active_radiation']
        """[W_{PAR} m-2ground] dictionary of incident photosynthetically active radiation.

        Notes:
            This input must be a dictionary having two keys, respectively 'sunlit' and 'shaded'
        """

        self.solar_inclination = self._inputs['solar_inclination']
        """(Rad) the angle between solar beam and the horizon.
        """
