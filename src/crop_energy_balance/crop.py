from pathlib import Path

from crop_energy_balance import utils
from crop_energy_balance.formalisms import canopy, weather, leaf, lumped_leaves, sunlit_shaded_leaves, component, soil
from crop_energy_balance.inputs import LumpedInputs, SunlitShadedInputs
from crop_energy_balance.params import Params, Constants

constants = Constants()


class CanopyStateVariables:
    def __init__(self,
                 inputs: LumpedInputs or SunlitShadedInputs):
        self.vapor_pressure_deficit = inputs.vapor_pressure_deficit
        self.vapor_pressure_slope = weather.calc_vapor_pressure_slope(
            utils.convert_kelvin_to_celsius(inputs.air_temperature, constants.absolute_zero))
        self.zero_displacement_height = canopy.calc_zero_displacement_height(
            canopy_height=inputs.canopy_height)
        self.roughness_length_for_momentum = canopy.calc_canopy_roughness_length_for_momentum(
            canopy_height=inputs.canopy_height)
        self.roughness_length_for_heat_transfer = canopy.calc_canopy_roughness_length_for_heat_transfer(
            canopy_height=inputs.canopy_height)
        self.wind_speed_at_canopy_height = canopy.calc_wind_speed_at_canopy_height(
            wind_speed=inputs.wind_speed,
            canopy_height=inputs.canopy_height,
            measurement_height=inputs.measurement_height)
        self.turbulent_diffusivity = canopy.calc_turbulent_diffusivity(
            von_karman_constant=constants.von_karman,
            wind_speed=inputs.wind_speed,
            canopy_height=inputs.canopy_height,
            zero_displacement_height=self.zero_displacement_height,
            canopy_roughness_length_for_momentum=self.roughness_length_for_heat_transfer,
            measurement_height=inputs.measurement_height)
        self.aerodynamic_resistance = canopy.calc_canopy_aerodynamic_resistance(
            wind_speed=inputs.wind_speed,
            canopy_height=inputs.canopy_height,
            reference_height=inputs.measurement_height,
            von_karman_constant=constants.von_karman)
        self.lumped_aerodynamic_resistance = canopy.calc_canopy_lumped_aerodynamic_resistance(
            canopy_aerodynamic_resistance=self.aerodynamic_resistance,
            vapor_pressure_slope=self.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

        self.sum_composed_conductances = None
        self.penman_energy = None
        self.total_penman_monteith_evaporative_energy = None
        self.source_temperature = None

    def calc_total_composed_conductances(self,
                                         crop_components: dict):
        self.sum_composed_conductances = sum(
            [1.0 / crop_component.composed_resistance for crop_component in crop_components.values()])

    def calc_total_evaporative_energy(self,
                                      crop_components: dict):
        self.penman_energy = canopy.calc_penman_evaporative_energy(
            canopy_aerodynamic_resistance=self.aerodynamic_resistance,
            canopy_net_radiation=sum([crop_component.net_radiation for crop_component in crop_components.values()]),
            vapor_pressure_slope=self.vapor_pressure_slope,
            vapor_pressure_deficit=self.vapor_pressure_deficit,
            psychrometric_constant=constants.psychrometric_constant,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity)
        self.total_penman_monteith_evaporative_energy = canopy.calc_penman_monteith_evaporative_energy(
            components_indices=list(crop_components.keys()),
            canopy_lumped_aerodynamic_resistance=self.lumped_aerodynamic_resistance,
            penman_evaporative_energy=self.penman_energy,
            composed_boundary_and_surface_conductances={component_index: crop_component.composed_conductance for
                                                        component_index, crop_component in crop_components.items()},
            net_radiation_fluxes={component_index: crop_component.net_radiation for
                                  component_index, crop_component in crop_components.items()},
            boundary_layer_resistances={component_index: crop_component.boundary_resistance for
                                        component_index, crop_component in crop_components.items()},
            vapor_pressure_slope=self.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

    def calc_source_temperature(self,
                                crop_components: dict,
                                inputs: LumpedInputs or SunlitShadedInputs):
        self.source_temperature = max(
            constants.absolute_zero,
            canopy.calc_temperature(
                air_temperature=inputs.air_temperature,
                canopy_aerodynamic_resistance=self.aerodynamic_resistance,
                canopy_net_radiation=sum([crop_component.net_radiation for crop_component in crop_components.values()]),
                penman_monteith_evaporative_energy=self.total_penman_monteith_evaporative_energy,
                air_density=constants.air_density,
                air_specific_heat_capacity=constants.air_specific_heat_capacity))


class Component:
    def __init__(self,
                 index: int):
        self.index = index

        self.absorbed_irradiance = None
        self.upper_cumulative_leaf_area_index = None
        self.lower_cumulative_leaf_area_index = None
        self.surface_resistance = None
        self.boundary_resistance = None
        self.composed_resistance = None
        self.composed_conductance = None
        self.net_longwave_radiation = None
        self.net_radiation = None
        self.penman_monteith_evaporative_energy = None
        self._temperature = None
        self.temperature = None

    def init_state_variables(self,
                             inputs: LumpedInputs or SunlitShadedInputs,
                             params: Params,
                             canopy_state_variables: CanopyStateVariables):
        self._temperature = inputs.air_temperature
        self.composed_resistance = component.calc_composed_resistance(
            surface_resistance=self.surface_resistance,
            boundary_layer_resistance=self.boundary_resistance,
            vapor_pressure_slope=canopy_state_variables.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.net_longwave_radiation = component.calc_net_longwave_radiation(
            air_temperature=inputs.air_temperature,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            atmospheric_emissivity=inputs.atmospheric_emissivity,
            extinction_coefficient=params.simulation.diffuse_black_extinction_coefficient,
            stefan_boltzman_constant=constants.stefan_boltzmann)
        self.net_radiation = component.calc_net_radiation(
            net_shortwave_radiation=utils.convert_photosynthetically_active_radition_into_global_radiation(
                self.absorbed_irradiance),
            net_longwave_radiation=self.net_longwave_radiation,
            is_soil=False)

    def calc_composed_conductance(self,
                                  canopy_state_variables: CanopyStateVariables):
        self.composed_conductance = component.calc_composed_conductance(
            composed_boundary_and_surface_resistance=self.composed_resistance,
            sum_composed_boundary_and_surface_conductances=canopy_state_variables.sum_composed_conductances,
            canopy_lumped_aerodynamic_resistance=canopy_state_variables.lumped_aerodynamic_resistance)

    def calc_evaporative_energy(self,
                                canopy_state_variables: CanopyStateVariables):
        self.penman_monteith_evaporative_energy = component.calc_evaporative_energy(
            net_radiation=self.net_radiation,
            boundary_layer_resistance=self.boundary_resistance,
            lumped_boundary_and_surface_resistance=self.composed_resistance,
            canopy_lumped_aerodynamic_resistance=canopy_state_variables.lumped_aerodynamic_resistance,
            penman_evaporative_energy=canopy_state_variables.penman_energy,
            penman_monteith_evaporative_energy=canopy_state_variables.total_penman_monteith_evaporative_energy,
            vapor_pressure_slope=canopy_state_variables.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

    def calc_temperature(self,
                         canopy_state_variables: CanopyStateVariables):
        self.temperature = component.calc_temperature(
            canopy_temperature=canopy_state_variables.source_temperature,
            boundary_layer_resistance=self.boundary_resistance,
            component_net_radiation=self.net_radiation,
            component_evaporative_energy=self.penman_monteith_evaporative_energy,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity)

    def update_temperature(self,
                           params: Params):
        previous_value = self._temperature
        self._temperature = self.temperature
        self.temperature += utils.calc_temperature_step(
            previous_value=previous_value,
            actual_value=self.temperature,
            step_fraction=params.numerical_resolution.step_fraction)


class SoilComponent(Component):
    def __init__(self,
                 index: int,
                 lower_cumulative_leaf_area_index: float):
        Component.__init__(self,
                           index=index)

        self.lower_cumulative_leaf_area_index = lower_cumulative_leaf_area_index

    def init_state_variables(self,
                             inputs: LumpedInputs or SunlitShadedInputs,
                             params: Params,
                             canopy_state_variables: CanopyStateVariables):
        self.absorbed_irradiance = inputs.absorbed_irradiance[self.index]['lumped']
        self.surface_resistance = soil.calc_surface_resistance(
            soil_saturation_ratio=inputs.soil_saturation_ratio,
            shape_parameter_1=params.simulation.soil_resistance_to_vapor_shape_parameter_1,
            shape_parameter_2=params.simulation.soil_resistance_to_vapor_shape_parameter_2)
        self.boundary_resistance = soil.calc_boundary_resistance(
            canopy_height=inputs.canopy_height,
            wind_speed=inputs.wind_speed,
            measurement_height=inputs.measurement_height,
            soil_roughness_length_for_momentum=params.simulation.soil_roughness_length_for_momentum,
            shape_parameter=params.simulation.soil_aerodynamic_resistance_shape_parameter,
            von_karman_constant=constants.von_karman)

        Component.init_state_variables(self,
                                       inputs=inputs,
                                       params=params,
                                       canopy_state_variables=canopy_state_variables)


class LeafComponent(Component):
    def __init__(self,
                 index: int,
                 upper_cumulative_leaf_area_index: float,
                 thickness: float):
        Component.__init__(self,
                           index=index)

        self.upper_cumulative_leaf_area_index = upper_cumulative_leaf_area_index
        self.lower_cumulative_leaf_area_index = upper_cumulative_leaf_area_index + thickness

        self.stomatal_sensibility = None


class LumpedLeafComponent(LeafComponent):
    def __init__(self, **kwargs):
        LeafComponent.__init__(self, **kwargs)

    def init_state_variables(self,
                             inputs: LumpedInputs,
                             params: Params,
                             canopy_state_variables: CanopyStateVariables):
        self.absorbed_irradiance = inputs.absorbed_irradiance[self.index]['lumped']
        self.stomatal_sensibility = leaf.calc_stomatal_sensibility(
            inputs.vapor_pressure_deficit,
            params.simulation.vapor_pressure_deficit_coefficient)
        self.surface_resistance = lumped_leaves.calc_leaf_layer_surface_resistance_to_vapor(
            absorbed_irradiance=self.absorbed_irradiance,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            stomatal_sensibility_to_water_status=self.stomatal_sensibility,
            global_extinction_coefficient=params.simulation.global_extinction_coefficient,
            maximum_stomatal_conductance=params.simulation.maximum_stomatal_conductance,
            residual_stomatal_conductance=params.simulation.residual_stomatal_conductance,
            shape_parameter=params.simulation.absorbed_par_50)
        self.boundary_resistance = lumped_leaves.calc_leaf_layer_boundary_resistance_to_vapor(
            wind_speed_at_canopy_height=inputs.wind_speed_at_canopy_height / 3600.0,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            wind_speed_extinction_coefficient=params.simulation.wind_speed_extinction_coef,
            characteristic_length=params.simulation.leaf_characteristic_length,
            shape_parameter=params.simulation.leaf_boundary_layer_shape_parameter,
            stomatal_density_factor=params.simulation.stomatal_density_factor)

        Component.init_state_variables(self,
                                       inputs=inputs,
                                       params=params,
                                       canopy_state_variables=canopy_state_variables)


class SunlitShadedLeafComponent(LeafComponent):
    def __init__(self, leaves_category: str, **kwargs):
        LeafComponent.__init__(self, **kwargs)

        self.leaves_category = leaves_category

    def init_state_variables(self,
                             inputs: SunlitShadedInputs,
                             params: Params,
                             canopy_state_variables: CanopyStateVariables):
        self.absorbed_irradiance = inputs.absorbed_irradiance[self.index][self.leaves_category]
        self.stomatal_sensibility = leaf.calc_stomatal_sensibility(
            inputs.vapor_pressure_deficit,
            params.simulation.vapor_pressure_deficit_coefficient)
        self.surface_resistance = sunlit_shaded_leaves.calc_leaf_layer_surface_resistance_to_vapor(
            leaves_category=self.leaves_category,
            incident_direct_irradiance=inputs.incident_par['direct'],
            incident_diffuse_irradiance=inputs.incident_par['diffuse'],
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            stomatal_sensibility_to_water_status=self.stomatal_sensibility,
            leaf_scattering_coefficient=params.simulation.leaf_scattering_coefficient,
            canopy_reflectance_to_direct_irradiance=params.simulation.canopy_reflectance_to_direct_irradiance,
            canopy_reflectance_to_diffuse_irradiance=params.simulation.canopy_reflectance_to_diffuse_irradiance,
            direct_extinction_coefficient=params.simulation.direct_extinction_coefficient,
            direct_black_extinction_coefficient=params.simulation.direct_black_extinction_coefficient,
            diffuse_extinction_coefficient=params.simulation.diffuse_extinction_coefficient,
            maximum_stomatal_conductance=params.simulation.maximum_stomatal_conductance,
            residual_stomatal_conductance=params.simulation.residual_stomatal_conductance,
            shape_parameter=params.simulation.absorbed_par_50,
            sublayers_number=params.simulation.sublayers_number)
        self.boundary_resistance = sunlit_shaded_leaves.calc_leaf_layer_boundary_resistance_to_vapor(
            leaves_category=self.leaves_category,
            wind_speed_at_canopy_height=inputs.wind_speed_at_canopy_height,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            direct_black_extinction_coefficient=params.simulation.direct_black_extinction_coefficient,
            wind_speed_extinction_coefficient=params.simulation.wind_speed_extinction_coef,
            characteristic_length=params.simulation.leaf_characteristic_length,
            shape_parameter=params.simulation.leaf_boundary_layer_shape_parameter,
            stomatal_density_factor=params.simulation.stomatal_density_factor)

        Component.init_state_variables(self,
                                       inputs=inputs,
                                       params=params,
                                       canopy_state_variables=canopy_state_variables)


class Canopy(dict):
    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __init__(self,
                 leaves_category: str,
                 inputs: LumpedInputs or SunlitShadedInputs = None,
                 params: Params = None,
                 inputs_path: Path = None,
                 params_path: Path = None):
        """Creates a class:`Shoot` object having either 'lumped' leaves or 'sunlit-shaded' leaves.
        Args:
            leaves_category: one of ('lumped', 'sunlit-shaded')
            inputs: see class`Inputs` and `SunlitShadedInputs`
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

        self.leaves_category = leaves_category
        if inputs_path:
            if self.leaves_category == 'lumped':
                self.inputs = LumpedInputs(inputs_path)
            else:
                self.inputs = SunlitShadedInputs(inputs_path)
        else:
            self.inputs = inputs

        if params_path:
            self.params = Params(params_path)
        else:
            self.params = params

        self.components_keys = self.inputs.components_keys

        self.state_variables = None

        self._set_components()

    def _set_components(self):
        """Sets canopy's components.
        """

        upper_cumulative_leaf_area_index = 0.0
        for index in reversed(self.components_keys):
            if index != -1:
                layer_thickness = self.inputs.leaf_layers[index]
                if self.leaves_category == 'lumped':
                    self[index] = LumpedLeafComponent(
                        index=index,
                        upper_cumulative_leaf_area_index=upper_cumulative_leaf_area_index,
                        thickness=layer_thickness)
                else:
                    self[index] = {
                        category: SunlitShadedLeafComponent(
                            leaves_category=category,
                            index=index,
                            upper_cumulative_leaf_area_index=upper_cumulative_leaf_area_index,
                            thickness=layer_thickness)
                        for category in ('sunlit', 'shaded')
                    }
                upper_cumulative_leaf_area_index += layer_thickness
            else:
                self[index] = SoilComponent(
                    index=index,
                    lower_cumulative_leaf_area_index=upper_cumulative_leaf_area_index)
