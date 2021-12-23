from pathlib import Path

from crop_energy_balance.crop import Crop, CropStateVariables
from crop_energy_balance.formalisms import canopy, weather
from crop_energy_balance.inputs import Inputs
from crop_energy_balance.params import Params, Constants
from crop_energy_balance.utils import is_almost_equal

constants = Constants()


class Solver:
    def __init__(self,
                 leaves_category: str,
                 inputs_path: Path = None,
                 inputs_dict: dict = None,
                 params_path: Path = None,
                 params_dict: dict = None):

        self.leaves_category = leaves_category

        self.inputs = Inputs(inputs_path=inputs_path) if inputs_path is not None else Inputs(inputs_dict=inputs_dict)
        self.params = Params(params_path=params_path) if params_path is not None else Params(params_dict=params_dict)
        self.params.update(inputs=self.inputs)

        self.crop = Crop(leaves_category=self.leaves_category, inputs=self.inputs, params=self.params)

        self.components = self.crop.extract_all_components()

        self.stability_iterations_number = 1
        self.iterations_number = 0
        self.error_temperature = None
        self.error_sensible_heat_flux = None

        self.init_state_variables()

        self.energy_balance = None

    def run(self, is_stability_considered=False):
        """Solves the steady-state energy balance.

        Args:
            is_stability_considered: If True then turbulence neutrality conditions are considered following
                Webber et al. (2016), otherwise False (default)

        References:
            Webber et al. (2016)
                Simulating canopy temperature for modelling heat stress in cereals.
                Environmental Modelling and Software 77, 143 - 155
        """
        # Solves the energy balance for neutral conditions
        self.solve_transient_energy_balance()

        # Corrects the energy balance for non-neutral conditions
        if is_stability_considered:
            is_acceptable_error = False
            while not is_acceptable_error and self.stability_iterations_number <= 100:
                self.stability_iterations_number += 1
                sensible_heat = self.crop.state_variables.sensible_heat_flux.copy()
                self.crop.state_variables.calc_aerodynamic_resistance(
                    inputs=self.crop.inputs, correct_stability=True)
                self.solve_transient_energy_balance()
                self.error_sensible_heat_flux = abs(self.crop.state_variables.sensible_heat_flux - sensible_heat)
                is_acceptable_error = is_almost_equal(actual=self.error_sensible_heat_flux, desired=0, decimal=2)

            if self.stability_iterations_number > 100:
                self.force_aerodynamic_resistance()
                self.solve_transient_energy_balance()

    def solve_transient_energy_balance(self):
        """Solves energy balance having fixed stability-related variables.
        """
        is_acceptable_error = False
        while not is_acceptable_error:
            self.iterations_number += 1
            self.update_state_variables()
            self.error_temperature = self.calc_error()
            self.update_temperature()
            self.calc_energy_balance()
            is_acceptable_error = self.determine_if_acceptable_error()

    def init_state_variables(self):
        self.crop.state_variables = CropStateVariables(inputs=self.crop.inputs, params=self.crop.params)
        for crop_component in self.components:
            crop_component.init_state_variables(self.crop.inputs, self.crop.params, self.crop.state_variables)

    def update_state_variables(self):
        self.crop.state_variables.calc_total_composed_conductances(crop_components=self.components)
        for crop_component in self.components:
            crop_component.calc_composed_conductance(self.crop.state_variables)

        self.crop.state_variables.calc_total_evaporative_energy(crop_components=self.components)
        for crop_component in self.components:
            crop_component.calc_evaporative_energy(self.crop.state_variables)

        self.crop.state_variables.calc_source_temperature(inputs=self.crop.inputs)

        self.crop.state_variables.calc_available_energy(crop_components=self.components)
        self.crop.state_variables.calc_net_radiation(soil_heat_flux=self.components[0].heat_flux)
        self.crop.state_variables.calc_sensible_heat_flux(inputs=self.crop.inputs)

        for crop_component in self.components:
            crop_component.calc_temperature(self.crop.state_variables)

    def update_temperature(self):
        for crop_component in self.components:
            crop_component.update_temperature(self.params)

    def calc_energy_balance(self):
        self.energy_balance = self.crop.state_variables.net_radiation - (
                self.crop.state_variables.total_penman_monteith_evaporative_energy +
                self.crop.state_variables.sensible_heat_flux +
                self.components[0].heat_flux)

    def force_aerodynamic_resistance(self):
        """Forces the aerodynamic resistance to a ratio of the neutral aerodynamic resistance, following
        Webber et al. (2016)

        References:
            Webber et al. (2016)
                Simulating canopy temperature for modelling heat stress in cereals.
                Environmental Modelling and Software 77, 143 - 155
        """
        neutral_friction_velocity = weather.calc_friction_velocity(
            wind_speed=self.crop.inputs.wind_speed,
            measurement_height=self.crop.inputs.measurement_height,
            zero_displacement_height=self.crop.state_variables.zero_displacement_height,
            roughness_length_for_momentum=self.crop.state_variables.roughness_length_for_momentum,
            stability_correction_for_momentum=0,
            von_karman_constant=constants.von_karman)
        neutral_aerodynamic_resistance = canopy.calc_aerodynamic_resistance(
            richardson_number=0,
            friction_velocity=neutral_friction_velocity,
            measurement_height=self.crop.inputs.measurement_height,
            zero_displacement_height=self.crop.state_variables.zero_displacement_height,
            roughness_length_for_heat=self.crop.state_variables.roughness_length_for_heat_transfer,
            stability_correction_for_heat=0,
            canopy_temperature=self.crop.state_variables.source_temperature,
            air_temperature=self.crop.inputs.air_temperature,
            von_karman_constant=constants.von_karman,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity)
        if self.crop.state_variables.sensible_heat_flux > 0:
            self.crop.state_variables.aerodynamic_resistance = 1.2 * neutral_aerodynamic_resistance
        else:
            self.crop.state_variables.aerodynamic_resistance = 0.8 * neutral_aerodynamic_resistance

    def calc_error(self) -> float:
        return sum(
            [abs(crop_component.temperature - crop_component._temperature) for crop_component in self.components])

    def determine_if_acceptable_error(self) -> bool:
        return self.error_temperature <= self.params.numerical_resolution.acceptable_temperature_error
