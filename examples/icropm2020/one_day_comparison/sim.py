from json import load
from math import radians
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from crop_irradiance.uniform_crops import (
    inputs as irradiance_inputs, params as irradiance_params, shoot as irradiance_canopy)

from crop_energy_balance import (
    inputs as eb_inputs, params as eb_params, crop as eb_canopy, solver as eb_solver)
from crop_energy_balance.examples.sources.demo import get_weather_data

with open('inputs.json', mode='r') as f:
    json_inputs = load(f)

with open('params.json', mode='r') as f:
    json_params = load(f)

leaf_layers = {3: 1.0,
               2: 1.0,
               1: 1.0,
               0: 1.0}


def get_irradiance_sim_inputs_and_params(
        is_bigleaf: bool,
        is_lumped: bool,
        incident_direct_par_irradiance: float,
        incident_diffuse_par_irradiance: float,
        solar_inclination_angle: float) -> (
        irradiance_inputs.LumpedInputs or irradiance_inputs.SunlitShadedInputs,
        irradiance_params.LumpedParams or irradiance_params.SunlitShadedInputs):
    vegetative_layers = {0: sum(leaf_layers.values())} if is_bigleaf else leaf_layers.copy()

    common_inputs = dict(
        leaf_layers=vegetative_layers,
        incident_direct_irradiance=incident_direct_par_irradiance,
        incident_diffuse_irradiance=incident_diffuse_par_irradiance,
        solar_inclination=radians(solar_inclination_angle))
    common_params = dict(
        leaf_reflectance=0.08,
        leaf_transmittance=0.07,
        leaves_to_sun_average_projection=0.5,
        sky_sectors_number=3,
        sky_type='soc',
        canopy_reflectance_to_diffuse_irradiance=0.057)
    if is_lumped:
        sim_inputs = irradiance_inputs.LumpedInputs(model='de_pury', **common_inputs)
        sim_params = irradiance_params.LumpedParams(model='de_pury', **common_params)
    else:
        sim_inputs = irradiance_inputs.SunlitShadedInputs(**common_inputs)
        sim_params = irradiance_params.SunlitShadedParams(**common_params)
    sim_params.update(sim_inputs)

    return sim_inputs, sim_params


def calc_absorbed_irradiance(
        is_bigleaf: bool,
        is_lumped: bool,
        incident_direct_par_irradiance: float,
        incident_diffuse_par_irradiance: float,
        solar_inclination_angle: float) -> dict:
    inputs, params = get_irradiance_sim_inputs_and_params(**locals())
    leaves_category = 'lumped' if is_lumped else 'sunlit-shaded'
    canopy = irradiance_canopy.Shoot(
        leaves_category=leaves_category,
        inputs=inputs,
        params=params)
    canopy.calc_absorbed_irradiance()

    absorbed_par_irradiance = {index: layer.absorbed_irradiance for index, layer in canopy.items()}

    absorbed_par_irradiance.update(
        {-1: {'lumped': sum([incident_direct_par_irradiance, incident_diffuse_par_irradiance]) - (
            sum([sum(v.absorbed_irradiance.values()) for v in canopy.values()]))}})

    return absorbed_par_irradiance


def get_energy_balance_inputs_and_params(
        vegetative_layers: dict,
        leaf_class_type: str,
        absorbed_par_irradiance: dict,
        actual_weather_data: pd.Series) -> (
        eb_inputs.LumpedInputs or eb_inputs.SunlitShadedInputs,
        eb_params.LumpedParams or eb_params.SunlitShadedParams):
    raw_inputs = json_inputs.copy()

    if leaf_class_type == 'lumped':
        incident_photosynthetically_active_radiation = {
            'lumped': sum(actual_weather_data[['incident_direct_irradiance', 'incident_diffuse_irradiance']])}
        solar_inclination_angle = None
    else:
        incident_photosynthetically_active_radiation = {
            'direct': actual_weather_data['incident_direct_irradiance'],
            'diffuse': actual_weather_data['incident_diffuse_irradiance']}
        solar_inclination_angle = actual_weather_data['solar_declination']

    raw_inputs.update(
        {"leaf_layers": vegetative_layers,
         "solar_inclination": solar_inclination_angle,
         "wind_speed": actual_weather_data['wind_speed'],
         "vapor_pressure": actual_weather_data['vapor_pressure'],
         "vapor_pressure_deficit": actual_weather_data['vapor_pressure_deficit'],
         "air_temperature": actual_weather_data['air_temperature'],
         "incident_photosynthetically_active_radiation": incident_photosynthetically_active_radiation,
         "absorbed_photosynthetically_active_radiation": absorbed_par_irradiance})

    if leaf_class_type == 'lumped':
        energy_balance_inputs = eb_inputs.LumpedInputs(raw_inputs)
        energy_balance_params = eb_params.LumpedParams(json_params)
    else:
        energy_balance_inputs = eb_inputs.SunlitShadedInputs(raw_inputs)
        energy_balance_params = eb_params.SunlitShadedParams(json_params)

    energy_balance_params.update(inputs=energy_balance_inputs)

    return energy_balance_inputs, energy_balance_params


def calc_temperature(
        vegetative_layers: dict,
        leaf_class_type: str,
        absorbed_par_irradiance: dict,
        actual_weather_data: pd.Series) -> dict:
    inputs, params = get_energy_balance_inputs_and_params(**locals())

    canopy = eb_canopy.Canopy(leaves_category=leaf_class_type, inputs=inputs, params=params)
    solver = eb_solver.Solver(canopy=canopy, inputs=inputs, params=params)
    solver.run()

    if leaf_class_type == 'lumped':
        return {index: {'lumped': layer.temperature} for index, layer in canopy.items()}
    else:
        res = {index: {'sunlit': layer['sunlit'].temperature, 'shaded': layer['shaded'].temperature}
               for index, layer in canopy.items() if index != -1}
        res.update({-1: {'lumped': canopy[-1].temperature}})
        return res


# def init_results_dict(leaves_category: str,
#                       canopy_layers: dict) -> dict:
#     value = {'lumped': []} if leaves_category == 'lumped' else {'sunlit': [], 'shaded': []}
#
#     return {k: value for k in canopy_layers.keys()}


def plot_irradiance_dynamic_comparison(incident_irradiance: pd.Series,
                                       all_cases_absorbed_irradiance: dict):
    cases = all_cases_absorbed_irradiance.keys()

    fig, axes = plt.subplots(ncols=len(cases), sharex='all', sharey='all', figsize=(15, 5))
    for i, case in enumerate(cases):
        plot_irradiance_dynamics(
            ax=axes[i],
            incident_par_irradiance=incident_irradiance,
            simulation_case=case,
            all_cases_data=all_cases_absorbed_irradiance)

    for i, ax in enumerate(axes):
        ax.grid()
        ax.legend()
        ax.set(xlabel='hour')
        if i == 0:
            ax.set_ylabel(r'$\mathregular{W_{PAR} \cdot m^{-2}_{ground}}$')

    fig.tight_layout()
    fig.savefig('figs/irradiance.png')
    plt.close()


def plot_temperature_dynamic_comparison(temperature_air: pd.Series,
                                        all_cases_temperature: dict):
    cases = all_cases_temperature.keys()

    fig, axes = plt.subplots(ncols=len(cases), sharex='all', sharey='all', figsize=(15, 5))
    for i, case in enumerate(cases):
        plot_temperature_dynamics(
            ax=axes[i],
            temperature_air=temperature_air,
            simulation_case=case,
            all_cases_data=all_cases_temperature)

    for i, ax in enumerate(axes):
        ax.legend()
        ax.set(xlabel='hour')
        if i == 0:
            ax.set_ylabel(r'$\mathregular{temperature\/[^\circ C]}$')

    fig.tight_layout()
    fig.savefig('figs/temperature.png')
    plt.close()


def get_summary_data(simulation_case: str,
                     all_cases_data: dict) -> dict:
    data_dynamic = all_cases_data[simulation_case]
    _, leaf_class = simulation_case.split('_')
    component_indexes = data_dynamic[0].keys()

    if leaf_class == 'lumped':
        summary_data = {k: [] for k in component_indexes}

        for hour in range(24):
            for component_index in component_indexes:
                summary_data[component_index].append(data_dynamic[hour][component_index]['lumped'])
    else:
        summary_data = {k: {'sunlit': [], 'shaded': []} for k in component_indexes if k != -1}
        summary_data.update({-1: []})

        for hour in range(24):
            for component_index in component_indexes:
                if component_index != -1:
                    summary_data[component_index]['sunlit'].append(
                        data_dynamic[hour][component_index]['sunlit'])
                    summary_data[component_index]['shaded'].append(
                        data_dynamic[hour][component_index]['shaded'])
                else:
                    summary_data[component_index].append(data_dynamic[hour][component_index]['lumped'])
    return summary_data


def plot_leaf_profile(vegetative_layers: {int, dict}):
    fig, axs = plt.subplots(ncols=4, figsize=(15, 5))
    for i, (k, v) in enumerate(vegetative_layers.items()):
        layer_indices = list(v.keys())
        axs[i].plot(list(v.values()), layer_indices, 'o-')
        axs[i].set(title=k, xlabel=r'$\mathregular{m^2_{leaf} \cdot m^{-2}_{ground}}$', ylabel='layer index [-]',
                   yticks=layer_indices)
    fig.tight_layout()
    fig.savefig('figs/layers.png')
    plt.close()


def plot_irradiance_dynamics(ax: plt.axis,
                             incident_par_irradiance: pd.Series,
                             simulation_case: str,
                             all_cases_data: dict):
    summary_data = get_summary_data(simulation_case, all_cases_data)
    _, leaf_class = simulation_case.split('_')
    component_indexes = summary_data.keys()

    ax.set_title(simulation_case.replace('_', ' '))
    ax.plot(range(24), incident_par_irradiance, label='incident', color='k', linestyle='--', linewidth=2)
    # ax.plot(range(24), abs_irradiance[-1], label=f'absorbed soil', color='brown', linewidth=2)

    if leaf_class == 'lumped':
        for component_index in component_indexes:
            if component_index != -1:
                ax.plot(range(24), summary_data[component_index], label=f'abs {component_index}')
    else:
        for component_index in component_indexes:
            if component_index != -1:
                ax.plot(range(24), summary_data[component_index]['sunlit'],
                        label=f'abs sunlit {component_index}')
                ax.plot(range(24), summary_data[component_index]['shaded'],
                        label=f'abs shaded {component_index}')

    pass


def plot_temperature_dynamics(ax: plt.axis,
                              temperature_air: pd.Series,
                              simulation_case: str,
                              all_cases_data: dict):
    summary_data = get_summary_data(simulation_case, all_cases_data)
    _, leaf_class = simulation_case.split('_')
    component_indexes = summary_data.keys()

    ax.set_title(simulation_case.replace('_', ' '))
    ax.plot(range(24), temperature_air, label='air', color='k', linestyle='--', linewidth=2)

    if leaf_class == 'lumped':
        for component_index in component_indexes:
            if component_index != -1:
                ax.plot(range(24), [(v - 273.15 if v > 273.15 else None) for v in summary_data[component_index]],
                        label=f'{component_index}')
    else:
        for component_index in component_indexes:
            if component_index != -1:
                ax.plot(range(24),
                        [(v - 273.15 if v > 273.15 else None) for v in summary_data[component_index]['sunlit']],
                        label=f'sunlit {component_index}')
                ax.plot(range(24),
                        [(v - 273.15 if v > 273.15 else None) for v in summary_data[component_index]['shaded']],
                        label=f'shaded {component_index}')

    pass


def plot_temperature_one_hour_comparison(hour: int,
                                         hourly_weather: pd.DataFrame,
                                         all_cases_absorbed_irradiance: dict,
                                         all_cases_temperature: dict):
    assert all_cases_temperature.keys() == all_cases_absorbed_irradiance.keys()

    cases = all_cases_temperature.keys()

    fig, axes = plt.subplots(nrows=2, ncols=len(cases), sharex='row', sharey='all', figsize=(15, 10))
    for i, case in enumerate(cases):
        plot_temperature_at_one_hour(
            ax=axes[1, i],
            hour=hour,
            temperature_air=hourly_weather.loc[:, 'air_temperature'],
            simulation_case=case,
            all_cases_data=all_cases_temperature)
        plot_irradiance_at_one_hour(
            ax=axes[0, i],
            hour=hour,
            incident_direct=hourly_weather.loc[:, 'incident_direct_irradiance'],
            incident_diffuse=hourly_weather.loc[:, 'incident_diffuse_irradiance'],
            simulation_case=case,
            all_cases_data=all_cases_absorbed_irradiance)

    for i, ax in enumerate(axes.flatten()):
        ax.legend()

    for ax in axes[:, 0]:
        ax.set_ylabel('Component index [-]')

    fig.tight_layout()
    fig.savefig('figs/temperature_at_one_hour.png')
    plt.close()


def plot_temperature_at_one_hour(ax: plt.axis,
                                 hour: int,
                                 temperature_air: pd.Series,
                                 simulation_case: str,
                                 all_cases_data: dict):
    summary_data = get_summary_data(simulation_case, all_cases_data)
    _, leaf_class = simulation_case.split('_')
    component_indexes = summary_data.keys()

    ax.axvline(temperature_air[hour], label='air', color='k', linestyle='--', linewidth=2)

    if leaf_class == 'lumped':
        y, x = zip(*[(i, summary_data[i][hour]) for i in component_indexes if i != -1])
        ax.plot([(v - 273.15 if v > 273.15 else None) for v in x], y, 'o-', label='lumped')
    else:
        y, x_sun, x_sh = zip(
            *[(i, summary_data[i]['sunlit'][hour], summary_data[i]['shaded'][hour]) for i in component_indexes if
              i != -1])

        ax.plot([(v - 273.15 if v > 273.15 else None) for v in x_sun], y, 'o-', color='y', label='sunlit')
        ax.plot([(v - 273.15 if v > 273.15 else None) for v in x_sh], y, 'o-', color='brown', label='shaded')

    ax.set_xlabel(r'$\mathregular{[^\circ C]}$')
    return


def plot_irradiance_at_one_hour(ax: plt.axis,
                                hour: int,
                                incident_direct: pd.Series,
                                incident_diffuse: pd.Series,
                                simulation_case: str,
                                all_cases_data: dict):
    summary_data = get_summary_data(simulation_case, all_cases_data)
    _, leaf_class = simulation_case.split('_')
    component_indexes = summary_data.keys()

    ax.set_title(simulation_case.replace('_', ' '))
    ax.axvline(incident_direct[hour], label='incident direct', color='y', linestyle='--', linewidth=2)
    ax.axvline(incident_diffuse[hour], label='incident diffuse', color='red', linestyle='--', linewidth=2)

    if leaf_class == 'lumped':
        y, x = zip(*[(i, summary_data[i][hour]) for i in component_indexes if i != -1])
        ax.plot(x, y, 'o-', label='lumped')
    else:
        y, x_sun, x_sh = zip(
            *[(i, summary_data[i]['sunlit'][hour], summary_data[i]['shaded'][hour]) for i in component_indexes if
              i != -1])

        ax.plot(x_sun, y, 'o-', color='y', label='sunlit')
        ax.plot(x_sh, y, 'o-', color='brown', label='shaded')

    ax.set_xlabel(r'$\mathregular{[W \cdot m^{-2}_{ground}]}$')
    return


if __name__ == '__main__':
    (Path(__file__).parent / 'figs').mkdir(exist_ok=True)
    weather_data = get_weather_data()

    irradiance = {}
    temperature = {}
    layers = {}

    for canopy_type, leaves_type in (('bigleaf', 'lumped'),
                                     ('bigleaf', 'sunlit-shaded'),
                                     ('layered', 'lumped'),
                                     ('layered', 'sunlit-shaded')):

        canopy_layers = {0: sum(leaf_layers.values())} if canopy_type == 'bigleaf' else leaf_layers
        layers.update({f'{canopy_type} {leaves_type}': canopy_layers})
        hourly_absorbed_irradiance = []
        hourly_temperature = []

        for date, w_data in weather_data.iterrows():
            incident_direct_irradiance = w_data['incident_direct_irradiance']
            incident_diffuse_irradiance = w_data['incident_diffuse_irradiance']
            solar_inclination = w_data['solar_declination']
            wind_speed = w_data['wind_speed']
            vapor_pressure_deficit = w_data['vapor_pressure_deficit']
            vapor_pressure = w_data['vapor_pressure']
            air_temperature = w_data['air_temperature']

            absorbed_irradiance = calc_absorbed_irradiance(
                is_bigleaf=(canopy_type == 'bigleaf'),
                is_lumped=(leaves_type == 'lumped'),
                incident_direct_par_irradiance=incident_direct_irradiance,
                incident_diffuse_par_irradiance=incident_diffuse_irradiance,
                solar_inclination_angle=solar_inclination)

            hourly_absorbed_irradiance.append(absorbed_irradiance)

            hourly_temperature.append(calc_temperature(
                vegetative_layers=canopy_layers,
                leaf_class_type=leaves_type,
                absorbed_par_irradiance=absorbed_irradiance,
                actual_weather_data=w_data))

        irradiance[f'{canopy_type}_{leaves_type}'] = hourly_absorbed_irradiance
        temperature[f'{canopy_type}_{leaves_type}'] = hourly_temperature

    plot_leaf_profile(vegetative_layers=layers)

    plot_irradiance_dynamic_comparison(
        incident_irradiance=(weather_data.loc[:, ['incident_direct_irradiance', 'incident_diffuse_irradiance']]).sum(
            axis=1),
        all_cases_absorbed_irradiance=irradiance)

    plot_temperature_dynamic_comparison(
        temperature_air=weather_data.loc[:, 'air_temperature'],
        all_cases_temperature=temperature)

    plot_temperature_one_hour_comparison(
        hour=12,
        hourly_weather=weather_data,
        all_cases_absorbed_irradiance=irradiance,
        all_cases_temperature=temperature)
