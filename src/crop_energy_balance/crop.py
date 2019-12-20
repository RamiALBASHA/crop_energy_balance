from crop_energy_balance import utils
from crop_energy_balance.formalisms import canopy, weather, leaf, lumped_leaves, component
from crop_energy_balance.inputs import LumpedInputs, SunlitShadedInputs
from crop_energy_balance.params import Params, Constants

constants = Constants()


class LeafLayer:
    def __init__(self,
                 index: int,
                 upper_cumulative_leaf_area_index: float,
                 thickness: float,
                 absorbed_irradiance: float):
        self.index = index
        self.upper_cumulative_leaf_area_index = upper_cumulative_leaf_area_index
        self.lower_cumulative_leaf_area_index = upper_cumulative_leaf_area_index + thickness
        self.absorbed_irradiance = absorbed_irradiance

        self.stomatal_sensibility = None
        self.surface_resistance = None
        self.boundary_resistance = None
        self.composed_resistance = None
        self.composed_conductance = None
        self.net_longwave_radiation = None
        self.net_radiation = None


class LumpedLeafLayer(LeafLayer):
    def __init__(self, **kwargs):
        LeafLayer.__init__(self, **kwargs)

    def init_state_variables(self,
                             inputs: LumpedInputs,
                             params: Params):
        self.stomatal_sensibility = leaf.calc_stomatal_sensibility(
            inputs.vapor_pressure_deficit,
            params.simulation.vapor_pressure_deficit_coefficient)
        self.surface_resistance = lumped_leaves.calc_leaf_layer_surface_resistance_to_vapor(
            absorbed_irradiance=inputs.incident_par,
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
        self.composed_resistance = component.calc_component_composed_resistance(
            surface_resistance=self.surface_resistance,
            boundary_layer_resistance=self.boundary_resistance,
            vapor_pressure_slope=inputs.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.net_longwave_radiation = component.calc_component_net_longwave_radiation(
            air_temperature=inputs.air_temperature,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            atmospheric_emissivity=inputs.atmospheric_emissivity,
            extinction_coefficient=params.simulation.diffuse_black_extinction_coefficient,
            stefan_boltzman_constant=constants.stefan_boltzmann)
        self.net_radiation = component.calc_component_net_radiation(
            net_shortwave_radiation=utils.convert_photosynthetically_active_radition_into_global_radiation(
                self.absorbed_irradiance),
            net_longwave_radiation=self.net_longwave_radiation,
            is_soil=False)

    def update_state_variables(self,
                               canopy_state_variable: CanopyStateVariables):
        self.composed_conductance = component.calc_component_composed_conductance(
            composed_boundary_and_surface_resistance=self.composed_resistance,
            sum_composed_boundary_and_surface_conductances=canopy_state_variable.sum_composed_conductance,
            canopy_lumped_aerodynamic_resistance=canopy_state_variable.lumped_aerodynamic_resistance)
        


class SunlitShadedLeafLayer(LeafLayer):
    def __init__(self, **kwargs):
        LeafLayer.__init__(self, **kwargs)


class CanopyStateVariables:
    def __init__(self,
                 inputs: LumpedInputs or SunlitShadedInputs):
        self.vapor_pressure_deficit = inputs.vapor_pressure_deficit

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
        self.vapor_pressure_slope = weather.calc_vapor_pressure_slope(
            utils.convert_kelvin_to_celsius(inputs.air_temperature, constants.absolute_zero))
        self.aerodynamic_resistance = canopy.calc_canopy_aerodynamic_resistance(
            wind_speed=inputs.wind_speed,
            canopy_height=inputs.canopy_height,
            reference_height=inputs.measurement_height,
            von_karman_constant=constants.von_karman)
        self.lumped_aerodynamic_resistance = canopy.calc_canopy_lumped_aerodynamic_resistance(
            canopy_aerodynamic_resistance=self.aerodynamic_resistance,
            vapor_pressure_slope=self.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

        self.penman_energy = None
        self.sum_composed_conductance = None
        self.total_penman_monteith_evaporative_energy = None

    def update(self,
               crop_components: dict):
        self.vapor_pressure_deficit = 1.  # TODO: replace by a function of source temperature
        self.sum_composed_conductance = sum(
            [1.0 / crop_component.composed_resistance for crop_component in crop_components.values()])
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


class Canopy(dict):
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

        self.state_variables = CanopyStateVariables(self.inputs)
        self.set_components(leaves_category)

        self.state_variables.update(
            inputs=self.inputs,
            params=self.params)
        self.calc_absorbed_irradiance()

    def set_components(self, leaves_category: str):
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

        pass

    def init_components_state_variables(self):
        for index in self._leaf_layer_indexes:
            self[index].init_state_variables(self.inputs, self.params)
