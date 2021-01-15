from crop_energy_balance.crop import Canopy, CanopyStateVariables
from crop_energy_balance.inputs import Inputs
from crop_energy_balance.params import Params
from crop_energy_balance.utils import is_almost_equal


class Solver:
    def __init__(self,
                 canopy: Canopy,
                 inputs: Inputs,
                 params: Params):
        self.canopy = canopy
        self.inputs = inputs
        self.params = params

        self.components = self.canopy.extract_all_components()

        self.stability_iterations_number = 1
        self.iterations_number = 0
        self.init_state_variables()

        self.energy_balance = None

    def run(self, correct_neutrality=False):
        self.solve_energy_balance()

        if correct_neutrality:
            is_acceptable_error = False
            while not is_acceptable_error:
                self.stability_iterations_number += 1
                if self.stability_iterations_number >= 100:
                    break
                sensible_heat = self.canopy.state_variables.sensible_heat_flux.copy()
                self.canopy.state_variables.update(inputs=self.inputs)
                self.solve_energy_balance()
                error = abs(self.canopy.state_variables.sensible_heat_flux - sensible_heat)
                is_acceptable_error = is_almost_equal(actual=error, desired=0, decimal=2)

    def solve_energy_balance(self):
        is_acceptable_error = False
        while not is_acceptable_error:
            self.iterations_number += 1
            self.update_state_variables()
            error = self.calc_error()
            self.update_temperature()
            self.calc_energy_balance()
            is_acceptable_error = self.determine_if_acceptable_error(error)

    def init_state_variables(self):
        self.canopy.state_variables = CanopyStateVariables(self.inputs)
        for crop_components in self.components:
            crop_components.init_state_variables(self.canopy.inputs, self.canopy.params, self.canopy.state_variables)

    def update_state_variables(self):
        self.canopy.state_variables.calc_total_composed_conductances(crop_components=self.components)
        for crop_component in self.components:
            crop_component.calc_composed_conductance(self.canopy.state_variables)

        self.canopy.state_variables.calc_total_evaporative_energy(crop_components=self.components)
        for crop_component in self.components:
            crop_component.calc_evaporative_energy(self.canopy.state_variables)

        self.canopy.state_variables.calc_source_temperature(inputs=self.canopy.inputs)

        self.canopy.state_variables.calc_available_energy(crop_components=self.components)
        self.canopy.state_variables.calc_net_radiation(soil_heat_flux=self.components[0].heat_flux)
        self.canopy.state_variables.calc_sensible_heat_flux(inputs=self.canopy.inputs)

        for crop_component in self.components:
            crop_component.calc_temperature(self.canopy.state_variables)

    def update_temperature(self):
        for crop_component in self.components:
            crop_component.update_temperature(self.params)

    def calc_energy_balance(self):
        self.energy_balance = self.canopy.state_variables.net_radiation - (
                self.canopy.state_variables.total_penman_monteith_evaporative_energy +
                self.canopy.state_variables.sensible_heat_flux +
                self.components[0].heat_flux)

    def calc_error(self) -> float:
        return sum(
            [abs(crop_component.temperature - crop_component._temperature) for crop_component in self.components])

    def determine_if_acceptable_error(self,
                                      error: float) -> bool:
        return error <= self.params.numerical_resolution.acceptable_temperature_error
