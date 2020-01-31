from json import load
from math import radians
import matplotlib.pyplot as plt

import pandas as pd
from crop_irradiance.uniform_crops import (
    inputs as irradiance_inputs, params as irradiance_params, shoot as irradiance_canopy)

from crop_energy_balance import (
    inputs as eb_inputs, params as eb_params, crop as eb_canopy, solver as eb_solver, utils)

with open('inputs.json', mode='r') as f:
    json_inputs = load(f)

with open('params.json', mode='r') as f:
    json_params = load(f)

leaf_layers = {3: 1.0,
               2: 1.0,
               1: 1.0,
               0: 1.0}


def get_irradiance_sim_inputs_and_params(
        leaf_scale_type: str,
        leaf_class_type: str,
        incident_direct_par_irradiance: float,
        incident_diffuse_par_irradiance: float,
        solar_inclination_angle: float) -> (
        irradiance_inputs.LumpedInputs or irradiance_inputs.SunlitShadedInputs,
        irradiance_params.LumpedParams or irradiance_params.SunlitShadedInputs):
    if leaf_scale_type == 'bigleaf':
        vegetative_layers = {0: sum(leaf_layers.values())}
    else:
        vegetative_layers = leaf_layers.copy()

    if leaf_class_type == 'lumped':
        sim_inputs = irradiance_inputs.LumpedInputs(
            leaf_layers=vegetative_layers,
            incident_irradiance=incident_direct_par_irradiance + incident_diffuse_par_irradiance)
        sim_params = irradiance_params.LumpedParams(
            extinction_coefficient=0.45)
    else:
        sim_inputs = irradiance_inputs.SunlitShadedInputs(
            leaf_layers=vegetative_layers,
            incident_direct_irradiance=incident_direct_par_irradiance,
            incident_diffuse_irradiance=incident_diffuse_par_irradiance,
            solar_inclination=radians(solar_inclination_angle))
        sim_params = irradiance_params.SunlitShadedParams(
            leaf_reflectance=0.07,
            leaf_transmittance=0.00,
            leaves_to_sun_average_projection=0.5,
            sky_sectors_number=3,
            sky_type='soc',
            canopy_reflectance_to_diffuse_irradiance=0.057)
        sim_params.update(sim_inputs)

    return sim_inputs, sim_params


def calc_absorbed_irradiance(
        leaf_scale_type: str,
        leaf_class_type: str,
        incident_direct_par_irradiance: float,
        incident_diffuse_par_irradiance: float,
        solar_inclination_angle: float) -> dict:
    inputs, params = get_irradiance_sim_inputs_and_params(**locals())
    canopy = irradiance_canopy.Shoot(
        leaves_category=leaf_class_type,
        inputs=inputs,
        params=params)

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
        inputs = eb_inputs.LumpedInputs(raw_inputs)
        params = eb_params.LumpedParams(json_params)
    else:
        inputs = eb_inputs.SunlitShadedInputs(raw_inputs)
        params = eb_params.SunlitShadedParams(json_params)

    params.update(inputs=inputs)

    return inputs, params


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


def get_weather_data(file_name: str) -> pd.DataFrame:
    raw_data = pd.read_csv(file_name, decimal='.', sep=';', skiprows=6).set_index('time')
    raw_data.loc[:, 'wind_speed'] = raw_data.apply(lambda x: x['wind_speed'] * 3600.0, axis=1)
    raw_data.loc[:, 'incident_direct_irradiance'] = raw_data['incident_global_irradiance'].apply(
        lambda x: utils.convert_global_irradiance_into_photosynthetically_active_radiation(x * 0.80))
    raw_data.loc[:, 'incident_diffuse_irradiance'] = raw_data['incident_global_irradiance'].apply(
        lambda x: utils.convert_global_irradiance_into_photosynthetically_active_radiation(x * 0.20))
    raw_data.loc[:, 'vapor_pressure_deficit'] = raw_data.apply(
        lambda x: utils.calc_vapor_pressure_deficit(x['air_temperature'], x['air_temperature'], x['relative_humidity']),
        axis=1)
    raw_data.loc[:, 'vapor_pressure'] = raw_data.apply(
        lambda x: x['vapor_pressure_deficit'] * x['relative_humidity'] / 100., axis=1)

    raw_data.drop(['incident_global_irradiance', 'relative_humidity'], axis=1, inplace=True)

    return raw_data


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
        ax.legend()
        ax.set(xlabel='hour')
        if i == 0:
            ax.set_ylabel(r'$\mathregular{W_{PAR} \cdot m^{-2}_{ground}}$')

    fig.tight_layout()
    fig.savefig('irradiance.png')
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
    fig.savefig('temperature.png')
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
                                         temperature_air: pd.Series,
                                         all_cases_temperature: dict):
    cases = all_cases_temperature.keys()

    fig, axes = plt.subplots(ncols=len(cases), sharex='all', sharey='all', figsize=(15, 5))
    for i, case in enumerate(cases):
        plot_temperature_at_one_hour(
            ax=axes[i],
            hour=hour,
            temperature_air=temperature_air,
            simulation_case=case,
            all_cases_data=all_cases_temperature)

    for i, ax in enumerate(axes):
        ax.legend()
        ax.set(xlabel='hour')
        if i == 0:
            ax.set_ylabel(r'$\mathregular{temperature\/[^\circ C]}$')

    fig.tight_layout()
    fig.savefig('temperature_at_one_hour.png')
    plt.close()


def plot_temperature_at_one_hour(ax: plt.axis,
                                 hour: int,
                                 temperature_air: pd.Series,
                                 simulation_case: str,
                                 all_cases_data: dict):
    summary_data = get_summary_data(simulation_case, all_cases_data)
    _, leaf_class = simulation_case.split('_')
    component_indexes = summary_data.keys()

    ax.set_title(simulation_case.replace('_', ' '))
    ax.axvline(temperature_air[hour], label='air', color='k', linestyle='--', linewidth=2)

    if leaf_class == 'lumped':
        y, x = zip(*[(i, summary_data[i][hour]) for i in component_indexes if i != -1])
        ax.plot([(v - 273.15 if v > 273.15 else None) for v in x], y, 'o-', label='lumped')
    else:
        y, x_sun, x_sh = zip(
            *[(i, summary_data[i]['sunlit'][hour], summary_data[i]['shaded'][hour]) for i in component_indexes if
              i != -1])

        ax.plot([(v - 273.15 if v > 273.15 else None) for v in x_sh], y, 'o-', label='shaded')
        ax.plot([(v - 273.15 if v > 273.15 else None) for v in x_sun], y, 'o-', label='sunlit')

    pass


def plot_weather(actual_wather: pd.DataFrame):
    fig, ((ax_irradiance, ax_temperature), (ax_wind_speed, ax_vpd)) = plt.subplots(ncols=2, nrows=2, sharex='all')

    day_hours = range(24)
    ax_irradiance.plot(
        day_hours, actual_wather.loc[:, ['incident_direct_irradiance', 'incident_diffuse_irradiance']].sum(axis=1))
    ax_irradiance.set_ylabel(r'$\mathregular{irradiance\/[W_{PAR} \cdot m^{-2}_{ground}]}$')

    ax_temperature.plot(day_hours, actual_wather.loc[:, 'air_temperature'])
    ax_temperature.set_ylabel(r'$\mathregular{temperature\/[^\circ C]}$')

    ax_wind_speed.plot(day_hours, actual_wather.loc[:, 'wind_speed'] / 3600.)
    ax_wind_speed.set_ylabel(r'$\mathregular{wind\/speed\/[m \cdot s^{-1}]}$')

    ax_vpd.plot(day_hours, actual_wather.loc[:, 'vapor_pressure_deficit'])
    ax_vpd.set_ylabel('VPD [kPa]')

    fig.tight_layout()
    fig.savefig('weather.png')
    plt.close()


if __name__ == '__main__':
    weather_data = get_weather_data('weather.csv')

    irradiance = {}
    temperature = {}

    for canopy_type, leaves_type in (('bigleaf', 'lumped'),
                                     ('bigleaf', 'sunlit-shaded'),
                                     ('layered', 'lumped'),
                                     ('layered', 'sunlit-shaded')):

        canopy_layers = {0: sum(leaf_layers.values())} if canopy_type == 'bigleaf' else leaf_layers
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
                leaf_scale_type=canopy_type,
                leaf_class_type=leaves_type,
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

    plot_weather(weather_data)

    plot_irradiance_dynamic_comparison(
        incident_irradiance=(weather_data.loc[:, ['incident_direct_irradiance', 'incident_diffuse_irradiance']]).sum(
            axis=1),
        all_cases_absorbed_irradiance=irradiance)

    plot_temperature_dynamic_comparison(
        temperature_air=weather_data.loc[:, 'air_temperature'],
        all_cases_temperature=temperature)

    plot_temperature_one_hour_comparison(
        hour=12,
        temperature_air=weather_data.loc[:, 'air_temperature'],
        all_cases_temperature=temperature)
