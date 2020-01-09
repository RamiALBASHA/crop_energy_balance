from crop_energy_balance.crop import Canopy, CanopyStateVariables
from crop_energy_balance.inputs import LumpedInputs, SunlitShadedInputs
from crop_energy_balance.params import LumpedParams, SunlitShadedParams


class Solver:
    def __init__(self,
                 canopy: Canopy,
                 inputs: LumpedInputs or SunlitShadedInputs,
                 params: LumpedParams or SunlitShadedParams):
        self.canopy = canopy
        self.inputs = inputs
        self.params = params

        self.components = self.canopy.extract_all_components()

        self.is_acceptable_error = False

        self.init_state_variables()

    def run(self):
        is_acceptable_error = False
        while not is_acceptable_error:
            self.update_state_variables()
            error = self.calc_error()
            self.update_temperature()
            is_acceptable_error = self.determine_if_acceptable_error(error)

    def init_state_variables(self):
        self.canopy.state_variables = CanopyStateVariables(self.inputs)
        for crop_components in self.components:
            crop_components.init_state_variables(self.canopy.inputs, self.canopy.params, self.canopy.state_variables)

    def update_state_variables(self):
        self.canopy.state_variables.calc_total_composed_conductances(crop_components=self.components)
        for crop_components in self.components:
            crop_components.calc_composed_conductance(self.canopy.state_variables)

        self.canopy.state_variables.calc_total_evaporative_energy(crop_components=self.components)
        for crop_components in self.components:
            crop_components.calc_evaporative_energy(self.canopy.state_variables)

        self.canopy.state_variables.calc_source_temperature(crop_components=self.components,
                                                            inputs=self.canopy.inputs)
        for crop_components in self.components:
            crop_components.calc_temperature(self.canopy.state_variables)

    def update_temperature(self):
        for crop_components in self.components:
            crop_components.update_temperature(self.params)

    def calc_error(self) -> float:
        return sum(
            [abs(crop_component.temperature - crop_component._temperature) for crop_component in self.components])

    def determine_if_acceptable_error(self,
                                      error: float) -> bool:
        return error <= self.params.numerical_resolution.acceptable_temperature_error
