from crop_energy_balance.crop import Canopy, CanopyStateVariables
from crop_energy_balance.inputs import Inputs
from crop_energy_balance.params import Params


class Solver:
    def __init__(self,
                 canopy: Canopy,
                 inputs: Inputs,
                 params: Params):
        self.canopy = canopy
        self.inputs = inputs
        self.params = params

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
        for index in self.canopy.components_keys:
            self.canopy[index].init_state_variables(self.canopy.inputs, self.canopy.params, self.canopy.state_variables)

    def update_state_variables(self):
        self.canopy.state_variables.calc_total_composed_conductances(crop_components=self.canopy)
        for index in self.canopy.components_keys:
            self.canopy[index].calc_composed_conductance(self.canopy.state_variables)

        self.canopy.state_variables.calc_total_evaporative_energy(crop_components=self.canopy)
        for index in self.canopy.components_keys:
            self.canopy[index].calc_evaporative_energy(self.canopy.state_variables)

        self.canopy.state_variables.calc_source_temperature(crop_components=self.canopy,
                                                            inputs=self.canopy.inputs)
        for index in self.canopy.components_keys:
            self.canopy[index].calc_temperature(self.canopy.state_variables)

    def update_temperature(self):
        for index in self.canopy.components_keys:
            self.canopy[index].update_temperature(self.params)

    def calc_error(self) -> float:
        return sum(
            [abs(crop_component.temperature - crop_component._temperature) for crop_component in self.canopy.values()])

    def determine_if_acceptable_error(self,
                                      error: float) -> bool:
        return error <= self.params.numerical_resolution.acceptable_temperature_error
