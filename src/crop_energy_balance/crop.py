from inspect import signature
from pathlib import Path

from crop_energy_balance import utils
from crop_energy_balance.formalisms import (canopy, weather, leaf, lumped_leaves, sunlit_shaded_leaves, component, soil,
                                            config)
from crop_energy_balance.inputs import Inputs
from crop_energy_balance.params import Params, Constants

constants = Constants()


class CropStateVariables:
    def __init__(self,
                 inputs: Inputs,
                 params: Params):
        self.sensible_heat_flux = config.PRECISION
        self.stability_correction_for_momentum = 0.0
        self.stability_correction_for_heat = 0.0

        self.source_temperature = inputs.air_temperature

        self.incident_irradiance = inputs.incident_irradiance
        self.vapor_pressure_deficit = inputs.vapor_pressure_deficit
        self.vapor_pressure_slope = weather.calc_vapor_pressure_slope(
            weather.convert_kelvin_to_celsius(inputs.air_temperature, constants.absolute_zero))
        self.zero_displacement_height = canopy.calc_zero_displacement_height(
            canopy_height=inputs.canopy_height,
            leaf_area_index=sum(inputs.leaf_layers.values()),
            drag_coefficient=params.simulation.drag_coefficient)
        self.roughness_length_for_momentum = canopy.calc_roughness_length_for_momentum_transfer(
            soil_roughness_length_for_momentum=params.simulation.soil_roughness_length_for_momentum,
            zero_plan_displacement_height=self.zero_displacement_height,
            canopy_height=inputs.canopy_height,
            total_leaf_area_index=sum(inputs.leaf_layers.values()),
            drag_coefficient=params.simulation.drag_coefficient)
        self.roughness_length_for_heat_transfer = canopy.calc_roughness_length_for_heat_transfer(
            roughness_length_for_momentum_transfer=self.roughness_length_for_momentum,
            ratio_heat_to_momentum_roughness_lengths=params.simulation.ratio_heat_to_momentum_canopy_roughness_lengths)
        self.wind_speed_at_canopy_height = canopy.calc_wind_speed_at_canopy_height(
            wind_speed=inputs.wind_speed,
            canopy_height=inputs.canopy_height,
            measurement_height=inputs.measurement_height,
            zero_displacement_height=self.zero_displacement_height,
            roughness_length_for_momentum=self.roughness_length_for_momentum)
        self.friction_velocity = canopy.calc_friction_velocity(
            wind_speed=inputs.wind_speed,
            measurement_height=inputs.measurement_height,
            zero_displacement_height=self.zero_displacement_height,
            roughness_length_for_momentum=self.roughness_length_for_momentum,
            stability_correction_for_momentum=self.stability_correction_for_momentum,
            von_karman_constant=constants.von_karman)
        self.monin_obukhov_length = canopy.calc_monin_obukhov_length(
            surface_temperature=self.source_temperature,
            sensible_heat_flux=self.sensible_heat_flux,
            friction_velocity=self.friction_velocity,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity,
            gravitational_acceleration=constants.gravitational_acceleration,
            von_karman_constant=constants.von_karman)
        self.richardson_number = canopy.calc_richardson_number(
            is_stable=self.sensible_heat_flux < 0,
            measurement_height=inputs.measurement_height,
            zero_displacement_height=self.zero_displacement_height,
            monin_obukhov_length=self.monin_obukhov_length)
        self.aerodynamic_resistance = None
        self.lumped_aerodynamic_resistance = None
        self.calc_aerodynamic_resistance(
            inputs=inputs,
            params=params)

        self.net_longwave_radiation = canopy.calc_net_longwave_radiation(
            air_temperature=inputs.air_temperature,
            air_vapor_pressure=inputs.air_vapor_pressure,
            atmospheric_emissivity=params.simulation.atmospheric_emissivity,
            stefan_boltzmann_constant=constants.stefan_boltzmann)

        # Canopy available energy and net radiation are not initialized to None so that to initialize the soil heat flux
        self.available_energy = weather.convert_photosynthetically_active_radiation_into_global_radiation(
            sum(self.incident_irradiance.values())) + self.net_longwave_radiation
        self.net_radiation = self.available_energy

        self.sum_composed_conductances = None
        self.penman_energy = None
        self.total_penman_monteith_evaporative_energy = None

    def calc_total_composed_conductances(self,
                                         crop_components: list):
        self.sum_composed_conductances = sum(
            [1.0 / crop_component.composed_resistance for crop_component in crop_components])

    def calc_total_evaporative_energy(self, crop_components: list):
        self.penman_energy = canopy.calc_penman_evaporative_energy(
            canopy_aerodynamic_resistance=self.aerodynamic_resistance,
            canopy_available_energy=self.available_energy,
            vapor_pressure_slope=self.vapor_pressure_slope,
            vapor_pressure_deficit=self.vapor_pressure_deficit,
            psychrometric_constant=constants.psychrometric_constant,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity)
        self.total_penman_monteith_evaporative_energy = canopy.calc_penman_monteith_evaporative_energy(
            canopy_lumped_aerodynamic_resistance=self.lumped_aerodynamic_resistance,
            penman_evaporative_energy=self.penman_energy,
            composed_boundary_and_surface_conductances=[crop_component.composed_conductance for
                                                        crop_component in crop_components],
            available_energy_fluxes=[crop_component.available_energy for crop_component in crop_components],
            boundary_layer_resistances=[crop_component.boundary_resistance for crop_component in crop_components],
            vapor_pressure_slope=self.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

    def calc_source_temperature(self, inputs: Inputs):
        self.source_temperature = max(
            constants.absolute_zero,
            canopy.calc_temperature(
                air_temperature=inputs.air_temperature,
                canopy_aerodynamic_resistance=self.aerodynamic_resistance,
                canopy_available_energy=self.available_energy,
                penman_monteith_evaporative_energy=self.total_penman_monteith_evaporative_energy,
                air_density=constants.air_density,
                air_specific_heat_capacity=constants.air_specific_heat_capacity))

    def calc_available_energy(self,
                              crop_components: list):
        self.available_energy = sum(crop_component.available_energy for crop_component in crop_components)

    def calc_net_radiation(self,
                           soil_heat_flux):
        self.net_radiation = self.available_energy + soil_heat_flux

    def calc_sensible_heat_flux(self, inputs: Inputs):
        self.sensible_heat_flux = canopy.calc_sensible_heat_flux(
            source_temperature=self.source_temperature,
            air_temperature=inputs.air_temperature,
            aerodynamic_resistance=self.aerodynamic_resistance,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity)

    def calc_aerodynamic_resistance(self,
                                    inputs: Inputs,
                                    params: Params):
        self.aerodynamic_resistance = canopy.calc_aerodynamic_resistance(
            richardson_number=self.richardson_number,
            friction_velocity=self.friction_velocity,
            measurement_height=inputs.measurement_height,
            zero_displacement_height=self.zero_displacement_height,
            roughness_length_for_heat=self.roughness_length_for_heat_transfer,
            stability_correction_for_heat=self.stability_correction_for_heat,
            canopy_temperature=self.source_temperature,
            air_temperature=inputs.air_temperature,
            richardon_threshold_free_convection=params.simulation.richardon_threshold_free_convection,
            von_karman_constant=constants.von_karman,
            air_density=constants.air_density,
            air_specific_heat_capacity=constants.air_specific_heat_capacity,
            free_convection_shape_parameter=params.simulation.free_convection_shape_parameter)

        self.lumped_aerodynamic_resistance = canopy.calc_lumped_aerodynamic_resistance(
            canopy_aerodynamic_resistance=self.aerodynamic_resistance,
            vapor_pressure_slope=self.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)


class Component:
    def __init__(self,
                 index: int):
        self.index = index

        self.surface_fraction = 1.0

        self.absorbed_irradiance = None
        self.upper_cumulative_leaf_area_index = None
        self.lower_cumulative_leaf_area_index = None
        self.vapor_pressure_deficit = None
        self.water_potential = None
        self.surface_resistance = None
        self.boundary_forced_convection_resistance = None
        self.boundary_free_convection_resistance = None
        self.boundary_resistance = None
        self.composed_resistance = None
        self.composed_conductance = None
        self.net_longwave_radiation = None
        self.available_energy = None
        self.penman_monteith_evaporative_energy = None
        self._temperature = None
        self.temperature = None

    def init_state_variables(self,
                             inputs: Inputs,
                             params: Params,
                             crop_state_variables: CropStateVariables):
        self._temperature = inputs.air_temperature
        self.temperature = inputs.air_temperature

    def calc_composed_resistance(self,
                                 params: Params,
                                 crop_state_variables: CropStateVariables):
        self.composed_resistance = component.calc_composed_resistance(
            surface_resistance=self.surface_resistance,
            boundary_layer_resistance=self.boundary_resistance,
            vapor_pressure_slope=crop_state_variables.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant,
            stomatal_density_factor=params.simulation.stomatal_density_factor)

    def calc_composed_conductance(self,
                                  crop_state_variables: CropStateVariables):
        self.composed_conductance = component.calc_composed_conductance(
            composed_boundary_and_surface_resistance=self.composed_resistance,
            sum_composed_boundary_and_surface_conductances=crop_state_variables.sum_composed_conductances,
            canopy_lumped_aerodynamic_resistance=crop_state_variables.lumped_aerodynamic_resistance)

    def calc_evaporative_energy(self,
                                crop_state_variables: CropStateVariables):
        self.penman_monteith_evaporative_energy = component.calc_evaporative_energy(
            available_energy=self.available_energy,
            boundary_layer_resistance=self.boundary_resistance,
            lumped_boundary_and_surface_resistance=self.composed_resistance,
            canopy_lumped_aerodynamic_resistance=crop_state_variables.lumped_aerodynamic_resistance,
            penman_evaporative_energy=crop_state_variables.penman_energy,
            penman_monteith_evaporative_energy=crop_state_variables.total_penman_monteith_evaporative_energy,
            vapor_pressure_slope=crop_state_variables.vapor_pressure_slope,
            psychrometric_constant=constants.psychrometric_constant)

    def calc_temperature(self,
                         crop_state_variables: CropStateVariables):
        self.temperature = component.calc_temperature(
            canopy_temperature=crop_state_variables.source_temperature,
            boundary_layer_resistance=self.boundary_resistance,
            available_energy=self.available_energy,
            evaporative_energy=self.penman_monteith_evaporative_energy,
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
        self.heat_flux = None

    def init_state_variables(self,
                             inputs: Inputs,
                             params: Params,
                             crop_state_variables: CropStateVariables):
        Component.init_state_variables(self,
                                       inputs=inputs,
                                       params=params,
                                       crop_state_variables=crop_state_variables)

        self.absorbed_irradiance = inputs.absorbed_irradiance[self.index]['lumped']
        self.surface_resistance = soil.calc_surface_resistance(
            soil_saturation_ratio=inputs.soil_saturation_ratio,
            shape_parameter_1=params.simulation.soil_resistance_to_vapor_shape_parameter_1,
            shape_parameter_2=params.simulation.soil_resistance_to_vapor_shape_parameter_2)
        self.boundary_resistance = soil.calc_forced_convection_resistance(
            wind_speed=inputs.wind_speed,
            canopy_height=inputs.canopy_height,
            measurement_height=inputs.measurement_height,
            zero_displacement_height=crop_state_variables.zero_displacement_height,
            canopy_roughness_length_for_momentum=crop_state_variables.roughness_length_for_momentum,
            soil_roughness_length_for_momentum=params.simulation.soil_roughness_length_for_momentum,
            shape_parameter=params.simulation.soil_aerodynamic_resistance_shape_parameter,
            von_karman_constant=constants.von_karman)
        self.net_longwave_radiation = soil.calc_net_longwave_radiation(
            canopy_top_net_longwave_radiation=crop_state_variables.net_longwave_radiation,
            canopy_leaf_area_index=sum(inputs.leaf_layers.values()),
            diffuse_black_extinction_coefficient=params.simulation.diffuse_black_extinction_coefficient)
        self.heat_flux = soil.calc_heat_flux(
            net_above_ground_radiation=crop_state_variables.net_radiation,
            is_diurnal=sum(inputs.incident_irradiance.values()) > 0)
        self.available_energy = component.calc_available_energy(
            net_shortwave_radiation=weather.convert_photosynthetically_active_radiation_into_global_radiation(
                self.absorbed_irradiance),
            net_longwave_radiation=self.net_longwave_radiation,
            soil_heat_flux=self.heat_flux)
        Component.calc_composed_resistance(self,
                                           params=params,
                                           crop_state_variables=crop_state_variables)


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

    def set_stomatal_sensibility_args(self, model_params: dict):
        model_args = model_params.copy()
        sensibility_funcs = {
            'leuning': leaf.calc_stomatal_sensibility_leuning,
            'tuzet': leaf.calc_stomatal_sensibility_tuzet,
            'misson': leaf.calc_stomatal_sensibility_misson}

        for name, param_args in model_params.items():
            args = signature(sensibility_funcs[name]).parameters.keys()
            model_args[name].update({k: getattr(self, k) for k in args if k not in param_args})

        return model_args

    def calc_boundary_resistance(self):
        forced_convection_resistance = self.boundary_forced_convection_resistance
        free_convection_resistance = self.boundary_free_convection_resistance
        if free_convection_resistance is None:
            boundary_resistance = forced_convection_resistance
        else:
            boundary_resistance = component.calc_boundary_layer_resistance(
                forced_convection_resistance=forced_convection_resistance,
                free_convection_resistance=free_convection_resistance)
        self.boundary_resistance = boundary_resistance

    def calc_boundary_and_composed_conductances(self, **kwargs):
        pass


class LumpedLeafComponent(LeafComponent):
    def __init__(self, **kwargs):
        LeafComponent.__init__(self, **kwargs)

    def init_state_variables(self,
                             inputs: Inputs,
                             params: Params,
                             crop_state_variables: CropStateVariables):
        Component.init_state_variables(self,
                                       inputs=inputs,
                                       params=params,
                                       crop_state_variables=crop_state_variables)

        self.vapor_pressure_deficit = inputs.vapor_pressure_deficit
        self.water_potential = inputs.soil_water_potential
        self.absorbed_irradiance = inputs.absorbed_irradiance[self.index]['lumped']
        self.stomatal_sensibility = leaf.calc_stomatal_sensibility(
            self.set_stomatal_sensibility_args(model_params=params.simulation.stomatal_sensibility))
        self.surface_resistance = lumped_leaves.calc_leaf_layer_surface_resistance_to_vapor(
            incident_direct_irradiance=inputs.incident_irradiance['direct'],
            incident_diffuse_irradiance=inputs.incident_irradiance['diffuse'],
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
            sublayers_number=params.simulation.sublayers_number,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.boundary_forced_convection_resistance = lumped_leaves.calc_leaf_layer_forced_convection_resistance(
            wind_speed_at_canopy_height=crop_state_variables.wind_speed_at_canopy_height / 3600.0,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            wind_speed_extinction_coefficient=params.simulation.wind_speed_extinction_coef,
            characteristic_length=params.simulation.leaf_characteristic_length,
            shape_parameter=params.simulation.leaf_boundary_layer_shape_parameter,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.net_longwave_radiation = lumped_leaves.calc_leaf_layer_net_longwave_radiation(
            canopy_top_net_longwave_radiation=crop_state_variables.net_longwave_radiation,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            diffuse_black_extinction_coefficient=params.simulation.diffuse_black_extinction_coefficient)
        self.available_energy = component.calc_available_energy(
            net_shortwave_radiation=weather.convert_photosynthetically_active_radiation_into_global_radiation(
                self.absorbed_irradiance),
            net_longwave_radiation=self.net_longwave_radiation,
            soil_heat_flux=0)
        self.calc_boundary_and_composed_conductances(
            inputs=inputs,
            params=params,
            crop_state_variables=crop_state_variables)

    def calc_boundary_and_composed_conductances(self,
                                                inputs: Inputs,
                                                params: Params,
                                                crop_state_variables: CropStateVariables):
        self.boundary_free_convection_resistance = lumped_leaves.calc_leaf_layer_free_convection_resistance(
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            layer_temperature=self.temperature,
            air_temperature=inputs.air_temperature,
            heat_molecular_diffusivity=constants.molecular_diffusivity_water_vapor,
            characteristic_length=params.simulation.leaf_characteristic_length,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.calc_boundary_resistance()
        Component.calc_composed_resistance(self,
                                           params=params,
                                           crop_state_variables=crop_state_variables)


class SunlitShadedLeafComponent(LeafComponent):
    def __init__(self, leaves_category: str, **kwargs):
        LeafComponent.__init__(self, **kwargs)

        self.leaves_category = leaves_category

    def init_state_variables(self,
                             inputs: Inputs,
                             params: Params,
                             crop_state_variables: CropStateVariables):
        Component.init_state_variables(self,
                                       inputs=inputs,
                                       params=params,
                                       crop_state_variables=crop_state_variables)

        self.vapor_pressure_deficit = inputs.vapor_pressure_deficit
        self.water_potential = inputs.soil_water_potential

        self.surface_fraction = sunlit_shaded_leaves.calc_leaf_fraction_per_leaf_layer(
            leaves_category=self.leaves_category,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            direct_black_extinction_coefficient=params.simulation.direct_black_extinction_coefficient)
        self.absorbed_irradiance = inputs.absorbed_irradiance[self.index][self.leaves_category]
        self.stomatal_sensibility = leaf.calc_stomatal_sensibility(
            self.set_stomatal_sensibility_args(model_params=params.simulation.stomatal_sensibility))
        self.surface_resistance = sunlit_shaded_leaves.calc_leaf_layer_surface_resistance_to_vapor(
            leaves_category=self.leaves_category,
            incident_direct_irradiance=inputs.incident_irradiance['direct'],
            incident_diffuse_irradiance=inputs.incident_irradiance['diffuse'],
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
            sublayers_number=params.simulation.sublayers_number,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.boundary_forced_convection_resistance = sunlit_shaded_leaves.calc_leaf_layer_forced_convection_resistance(
            leaves_category=self.leaves_category,
            wind_speed_at_canopy_height=crop_state_variables.wind_speed_at_canopy_height / 3600.0,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            direct_black_extinction_coefficient=params.simulation.direct_black_extinction_coefficient,
            wind_speed_extinction_coefficient=params.simulation.wind_speed_extinction_coef,
            characteristic_length=params.simulation.leaf_characteristic_length,
            shape_parameter=params.simulation.leaf_boundary_layer_shape_parameter,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.net_longwave_radiation = sunlit_shaded_leaves.calc_leaf_layer_net_longwave_radiation(
            leaves_category=self.leaves_category,
            canopy_top_net_longwave_radiation=crop_state_variables.net_longwave_radiation,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            direct_black_extinction_coefficient=params.simulation.direct_black_extinction_coefficient,
            diffuse_black_extinction_coefficient=params.simulation.diffuse_black_extinction_coefficient)
        self.available_energy = component.calc_available_energy(
            net_shortwave_radiation=weather.convert_photosynthetically_active_radiation_into_global_radiation(
                self.absorbed_irradiance),
            net_longwave_radiation=self.net_longwave_radiation,
            soil_heat_flux=0)
        self.calc_boundary_and_composed_conductances(
            inputs=inputs,
            params=params,
            crop_state_variables=crop_state_variables)

    def calc_boundary_and_composed_conductances(self,
                                                inputs: Inputs,
                                                params: Params,
                                                crop_state_variables: CropStateVariables):
        self.boundary_free_convection_resistance = sunlit_shaded_leaves.calc_leaf_layer_free_convection_resistance(
            leaves_category=self.leaves_category,
            upper_cumulative_leaf_area_index=self.upper_cumulative_leaf_area_index,
            lower_cumulative_leaf_area_index=self.lower_cumulative_leaf_area_index,
            layer_temperature=self.temperature,
            air_temperature=inputs.air_temperature,
            heat_molecular_diffusivity=constants.molecular_diffusivity_water_vapor,
            direct_black_extinction_coefficient=params.simulation.direct_black_extinction_coefficient,
            characteristic_length=params.simulation.leaf_characteristic_length,
            stomatal_density_factor=params.simulation.stomatal_density_factor)
        self.calc_boundary_resistance()
        Component.calc_composed_resistance(self,
                                           params=params,
                                           crop_state_variables=crop_state_variables)


class Crop(dict):
    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __init__(self,
                 leaves_category: str,
                 inputs: Inputs = None,
                 params: Params = None,
                 inputs_path: Path = None,
                 params_path: Path = None):
        """Creates a class:`Canopy` object having either 'lumped' leaves or 'sunlit-shaded' leaves.
        Args:
            leaves_category: one of ('lumped', 'sunlit-shaded')
            inputs: see class`Inputs`
            params: see class`Params`
        Notes:
            The created canopy can implicitly be 'big-leaf' or a 'layered'. If the attribute `leaf_layers` of the
                :Class:`inputs` object has only one layer, then the resulting canopy is a 'big-leaf', otherwise if the
                dictionary has more than one key, then the canopy is 'layered'.
            Leaf layers indexes in `leaf_layers` must be ordered so that the youngest leaf layer has the highest index
                value, and inversely, the oldest leaf layer has the least value. Not respecting this order will
                definitely lead to erroneous calculations.
        """

        dict.__init__(self)

        self.leaves_category = leaves_category
        if inputs_path:
            self.inputs = Inputs(inputs_path=inputs_path)
        else:
            self.inputs = inputs

        if params_path:
            self.params = Params(params_path=params_path)
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

    def extract_all_components(self, ignore_nocturnal_sunlit: bool):
        """Returns a list of all canopy components."""
        if self.leaves_category == 'lumped':
            return [self[key] for key in self.components_keys]
        else:
            res = []
            for component_key in self.components_keys:
                if component_key != -1:
                    if not (ignore_nocturnal_sunlit and self.inputs.incident_irradiance['direct'] == 0):
                        res.append(self[component_key]['sunlit'])
                    res.append(self[component_key]['shaded'])
                else:
                    res.append(self[component_key])
            return res
