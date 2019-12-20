from crop_energy_balance.formalisms import canopy, weather, leaf, lumped_leaves, component
from crop_energy_balance.inputs import LumpedInputs, SunlitShadedInputs
from crop_energy_balance.params import Params, Constants

from crop_energy_balance import utils

constants = Constants()


class LeafLayer:
    def __init__(self,
                 index: int,
                 upper_cumulative_leaf_area_index: float,
                 thickness: float):
        self.index = index
        self.upper_cumulative_leaf_area_index = upper_cumulative_leaf_area_index
        self.lower_cumulative_leaf_area_index = upper_cumulative_leaf_area_index + thickness

        self.stomatal_closure_fraction = None
        self.surface_resistance = None
        self.boundary_resistance = None
        self.composed_resistance = None

    def calc_state_variables(self,
                             inputs: LumpedInputs or SunlitShadedInputs,
                             params: Params):
        self.stomatal_closure_fraction = leaf.calc_stomatal_sensibility(
            inputs.vapor_pressure_deficit,
            params.simulation.vapor_pressure_deficit_coefficient)


class LumpedLeafLayer(LeafLayer):
    def __init__(self, **kwargs):
        LeafLayer.__init__(self, **kwargs)

    def calc_state_variables(self,
                             inputs: LumpedInputs,
                             params: Params):
        self.stomatal_stomatal_sensibility = leaf.calc_stomatal_sensibility(
            inputs.vapor_pressure_deficit,
            params.simulation.vapor_pressure_deficit_coefficient)
        self.surface_resistance = lumped_leaves.calc_leaf_layer_surface_resistance_to_vapor(
            absorbed_irradiance=inputs.incident_par,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            stomatal_sensibility_to_water_status=self.stomatal_stomatal_sensibility,
            global_extinction_coefficient=params.simulation.global_extinction_coefficient,
            maximum_stomatal_conductance=params.simulation.maximum_stomatal_conductance,
            residual_stomatal_conductance=params.simulation.residual_stomatal_conductance,
            shape_parameter=params.simulation.absorbed_par_50)
        self.boundary_resistances = lumped_leaves.calc_leaf_layer_boundary_resistance_to_vapor(
            wind_speed_at_canopy_height=inputs.wind_speed_at_canopy_height / 3600.0,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            wind_speed_extinction_coefficient=params.simulation.wind_speed_extinction_coef,
            characteristic_length=params.simulation.leaf_characteristic_length,
            shape_parameter=params.simulation.leaf_boundary_layer_shape_parameter,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.composed_resistances = component.calc_component_composed_resistance(
            surface_resistance=self.surface_resistance,
            boundary_layer_resistance=self.boundary_resistance,
            vapor_pressure_slope=inputs.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant,
            stomatal_density_factor=params.simulation.stomatal_density_factor)


class SunlitShadedLeafLayer(LeafLayer):
    def __init__(self, **kwargs):
        LeafLayer.__init__(self, **kwargs)


class Shoot(dict):
    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __init__(self,
                 leaves_category: str,
                 inputs: LumpedInputs or SunlitShadedInputs,
                 params: Params):
        """Creates a class:`Shoot` object having either 'lumped' leaves or 'sunlit-shaded' leaves.
        Args:
            leaves_category: one of ('lumped', 'sunlit-shaded')
            inputs: see class`LumpedInputs` and `SunlitShadedInputs`
            params: see class`Params`
        Notes:
            The created shoot can implicitly be 'big-leaf' or a 'layered'. If the attribute `leaf_layers` of the
                :Class:`inputs` object has only one layer, then the resulting shoot is a 'big-leaf', otherwise if the
                dictionary has more than one key, then the shoot is 'layered'.
            Leaf layers indexes in `leaf_layers` must be ordered so that the youngest leaf layer has the highest index
                value, and inversely, the oldest leaf layer has the least value. Not respecting this order will
                definitely lead to erroneous calculations.
        """

        dict.__init__(self)

        self.inputs = inputs
        self.params = params
        self._leaf_layer_indexes = list(reversed(sorted(inputs.leaf_layers.keys())))

        self.zero_displacement_height = None
        self.canopy_roughness_length_for_momentum = None
        self.canopy_roughness_length_for_heat_transfer = None
        self.wind_speed_at_canopy_height = None
        self.turbulent_diffusivity = None
        self.canopy_aerodynamic_resistance = None

        self.set_leaf_layers(leaves_category)
        self.calc_absorbed_irradiance()

    def set_leaf_layers(self, leaves_category: str):
        """Sets leaf layers of the shoot.
        Args:
            leaves_category: one of ('lumped', 'sunlit-shaded')
        """

        upper_cumulative_leaf_area_index = 0.0
        for index in self._leaf_layer_indexes:
            layer_thickness = self.inputs.leaf_layers[index]
            if leaves_category == 'lumped':
                self[index] = LumpedLeafLayer(index=index,
                                              upper_cumulative_leaf_area_index=upper_cumulative_leaf_area_index,
                                              thickness=layer_thickness)
            else:
                self[index] = SunlitShadedLeafLayer(index=index,
                                                    upper_cumulative_leaf_area_index=upper_cumulative_leaf_area_index,
                                                    thickness=layer_thickness)

            upper_cumulative_leaf_area_index += layer_thickness

    def calc_absorbed_irradiance(self):
        """Calculates the absorbed irradiance by shoot's layers.
        """
        for index in self._leaf_layer_indexes:
            self[index].calc_absorbed_irradiance(self.inputs, self.params)

    def calc_canopy_state_variables(self):
        self.zero_displacement_height = canopy.calc_zero_displacement_height(
            canopy_height=self.inputs.canopy_height)
        self.canopy_roughness_length_for_momentum = canopy.calc_canopy_roughness_length_for_momentum(
            canopy_height=self.inputs.canopy_height)
        self.canopy_roughness_length_for_heat_transfer = canopy.calc_canopy_roughness_length_for_heat_transfer(
            canopy_height=self.inputs.canopy_height)
        self.wind_speed_at_canopy_height = canopy.calc_wind_speed_at_canopy_height(
            wind_speed=self.inputs.wind_speed,
            canopy_height=self.inputs.canopy_height,
            measurement_height=self.inputs.measurement_height)
        self.turbulent_diffusivity = canopy.calc_turbulent_diffusivity(
            von_karman_constant=constants.von_karman,
            wind_speed=self.inputs.wind_speed,
            canopy_height=self.inputs.canopy_height,
            zero_displacement_height=self.zero_displacement_height,
            canopy_roughness_length_for_momentum=self.canopy_roughness_length_for_heat_transfer,
            measurement_height=self.inputs.measurement_height)
        self.canopy_aerodynamic_resistance = canopy.calc_canopy_aerodynamic_resistance(
            wind_speed=self.inputs.wind_speed,
            canopy_height=self.inputs.canopy_height,
            reference_height=self.inputs.measurement_height,
            von_karman_constant=constants.von_karman)
        self.vapor_pressure_slope = weather.calc_vapor_pressure_slope(
            utils.convert_kelvin_to_celsius(self.inputs.air_temperature, constants.absolute_zero))
        self.canopy_lumped_aerodynamic_resistance = canopy.calc_canopy_lumped_aerodynamic_resistance(
            canopy_aerodynamic_resistance=self.canopy_aerodynamic_resistance,
            vapor_pressure_slope=self.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

        pass

    def calc_components_state_variables(self):
        pass
