from json import load
from pathlib import Path

from crop_energy_balance import params
from crop_energy_balance.formalisms import weather

constants = params.Constants()


class Inputs:
    def __init__(self,
                 inputs_dict: dict = None,
                 inputs_path: Path = None):
        if inputs_dict is None:
            with open(str(inputs_path), mode='r', encoding='utf-8') as f:
                inputs_dict = load(f)

        self._inputs = self._fmt_inputs(inputs_dict)

        self.measurement_height = self._inputs['measurement_height']
        """[m] height at which the meteorological variables are measured"""

        self.canopy_height = self._inputs['canopy_height']
        """[m] height of the canopy"""

        self.soil_saturation_ratio = self._inputs['soil_saturation_ratio']
        """[-] ratio of actual to potential volumetric water content in the soil"""

        self.air_temperature = weather.convert_celsius_to_kelvin(self._inputs['air_temperature'],
                                                                 constants.absolute_zero)
        """[K] temperature of the air"""

        self.wind_speed = max(0.01 * 3600., self._inputs['wind_speed'])
        """[m h-1] wind speed"""

        self.air_relative_humidity = self._inputs['relative_humidity']
        """[-] air relative humidity (%, between 0 and 100)"""

        self.air_vapor_pressure = self._inputs['vapor_pressure']
        """[kPa] vapor pressure of the air"""

        self.vapor_pressure_deficit = self._inputs['vapor_pressure_deficit']
        """[kPa] vapor pressure deficit of the air"""

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

        self.incident_irradiance = self._inputs['incident_photosynthetically_active_radiation']
        """[W_{PAR} m-2ground] dictionary of incident photosynthetically active radiation.

        Notes:
            This dictionary must provide direct and diffuse PAR with 'direct' and 'diffuse' keys, respectively.
        """

        self.absorbed_irradiance = self._inputs['absorbed_photosynthetically_active_radiation']
        """[W_{PAR} m-2ground] dictionary of absorbed photosynthetically active radiation per crop component.

        Notes:
            For lumped leaves, the absorbed irradiance per leaf layer is a dictionary having the key 'lumped'.
            For sunlit and shaded leaves, the absorbed irradiance per leaf layer is a dictionary having the keys
                'sunlit' and 'shaded'.
            The uppermost component must have the highest number while the lowermost component, i.e. soil, has the
                lowest number which must be equal to -1.
        """

        self.soil_water_potential = self._inputs[
            'soil_water_potential'] if 'soil_water_potential' in self._inputs else None
        """[MPa] soil water potential.
        """

        self.solar_inclination = self._inputs['solar_inclination']
        """(Rad) the angle between solar beam and the horizon.
        """

        self.components_keys = sorted(list(self.leaf_layers.keys()) + [-1])

    @staticmethod
    def _fmt_inputs(inputs: dict):
        for k in ('leaf_layers', 'absorbed_photosynthetically_active_radiation'):
            inputs[k] = {int(key): value for key, value in inputs[k].items()}
        return inputs
