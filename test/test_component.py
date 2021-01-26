from crop_energy_balance.formalisms import component
from crop_energy_balance.params import Constants
from crop_energy_balance.utils import assert_trend

constants = Constants()


def test_calc_composed_resistance():
    assert 0 == component.calc_composed_resistance(surface_resistance=0,
                                                   boundary_layer_resistance=0,
                                                   vapor_pressure_slope=1.0,
                                                   psychrometric_constant=constants.psychrometric_constant,
                                                   stomatal_density_factor=1)

    assert 1 == component.calc_composed_resistance(surface_resistance=1,
                                                   boundary_layer_resistance=0,
                                                   vapor_pressure_slope=1.0,
                                                   psychrometric_constant=constants.psychrometric_constant,
                                                   stomatal_density_factor=1)

    assert 2 == component.calc_composed_resistance(surface_resistance=0,
                                                   boundary_layer_resistance=1,
                                                   vapor_pressure_slope=constants.psychrometric_constant,
                                                   psychrometric_constant=constants.psychrometric_constant,
                                                   stomatal_density_factor=1)


def test_calc_composed_conductance():
    assert 1 == component.calc_composed_conductance(composed_boundary_and_surface_resistance=1,
                                                    sum_composed_boundary_and_surface_conductances=1,
                                                    canopy_lumped_aerodynamic_resistance=0)

    assert 0.5 == component.calc_composed_conductance(composed_boundary_and_surface_resistance=1,
                                                      sum_composed_boundary_and_surface_conductances=1,
                                                      canopy_lumped_aerodynamic_resistance=1)


def test_calc_evaporative_energy():
    assert 0 == component.calc_evaporative_energy(available_energy=0,
                                                  boundary_layer_resistance=0,
                                                  lumped_boundary_and_surface_resistance=1,
                                                  canopy_lumped_aerodynamic_resistance=0,
                                                  penman_evaporative_energy=1,
                                                  penman_monteith_evaporative_energy=1,
                                                  vapor_pressure_slope=1,
                                                  psychrometric_constant=constants.psychrometric_constant)

    assert 0 == component.calc_evaporative_energy(available_energy=0,
                                                  boundary_layer_resistance=0,
                                                  lumped_boundary_and_surface_resistance=1,
                                                  canopy_lumped_aerodynamic_resistance=1,
                                                  penman_evaporative_energy=1,
                                                  penman_monteith_evaporative_energy=1,
                                                  vapor_pressure_slope=1,
                                                  psychrometric_constant=constants.psychrometric_constant)

    assert 1 == component.calc_evaporative_energy(available_energy=1,
                                                  boundary_layer_resistance=1,
                                                  lumped_boundary_and_surface_resistance=1,
                                                  canopy_lumped_aerodynamic_resistance=0,
                                                  penman_evaporative_energy=1,
                                                  penman_monteith_evaporative_energy=1,
                                                  vapor_pressure_slope=constants.psychrometric_constant,
                                                  psychrometric_constant=constants.psychrometric_constant)


def test_calc_temperature():
    temperature = 25 + 273
    assert temperature == component.calc_temperature(canopy_temperature=temperature,
                                                     boundary_layer_resistance=1,
                                                     available_energy=1,
                                                     evaporative_energy=1,
                                                     air_density=constants.air_density,
                                                     air_specific_heat_capacity=constants.air_specific_heat_capacity)

    assert_trend(expected_trend='+',
                 values=[component.calc_temperature(canopy_temperature=temperature,
                                                    boundary_layer_resistance=ra,
                                                    available_energy=2,
                                                    evaporative_energy=1,
                                                    air_density=constants.air_density,
                                                    air_specific_heat_capacity=constants.air_specific_heat_capacity)
                         for ra in range(4)])

    assert_trend(expected_trend='+',
                 values=[component.calc_temperature(canopy_temperature=temperature,
                                                    boundary_layer_resistance=1,
                                                    available_energy=a,
                                                    evaporative_energy=1,
                                                    air_density=constants.air_density,
                                                    air_specific_heat_capacity=constants.air_specific_heat_capacity)
                         for a in range(4)])
