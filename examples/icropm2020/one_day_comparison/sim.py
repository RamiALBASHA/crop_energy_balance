import importlib.util
from json import load
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from crop_irradiance.uniform_crops import (
    inputs as irradiance_inputs, params as irradiance_params, shoot as irradiance_canopy)

from crop_energy_balance import (
    inputs as eb_inputs, params as eb_params, crop as eb_canopy, solver as eb_solver)

# import demo.py (outside src)
spec = importlib.util.spec_from_file_location(
    name="demo",
    location=r"../../../examples/sources/demo.py")

# spec = importlib.util.spec_from_file_location(
#     name="demo",
#     location=r"C:\Users\ralbasha\Documents\dvlp\crop_energy_balance_2\crop_energy_balance\examples\sources\demo.py")
foo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(foo)
foo.get_weather_data()
# from crop_energy_balance.examples.sources.demo import get_weather_data

with open('inputs.json', mode='r') as f:
    json_inputs = load(f)

with open('params.json', mode='r') as f:
    json_params = load(f)

leaf_layers = {3: 1.0,
               2: 1.0,
               1: 1.0,
               0: 1.0}

UNITS_MAP = {
    'net_radiation': (r'$\mathregular{R_n}$', r'$\mathregular{[W\/m^{-2}_{ground}]}$'),
    'sensible_heat_flux': ('H', r'$\mathregular{[W\/m^{-2}_{ground}]}$'),
    'total_penman_monteith_evaporative_energy': (r'$\mathregular{\lambda E}$', r'$\mathregular{[W\/m^{-2}_{ground}]}$'),
    'soil_heat_flux': ('G', r'$\mathregular{[W\/m^{-2}_{ground}]}$'),
    'energy_balance': ('balance', r'$\mathregular{[W\/m^{-2}_{ground}]}$'),
    'richardson_number': ('Ri', '[-]'),
    'stability_correction_for_momentum': (r'$\mathregular{\Phi_m}$', '[m]'),
    'stability_correction_for_heat': (r'$\mathregular{\Phi_h}$', '[m]'),
    'monin_obukhov_length': ('L', '[m]')
}


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
        solar_inclination=solar_inclination_angle)
    common_params = dict(
        leaf_reflectance=0.08,
        leaf_transmittance=0.07,
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
        absorbed_par_irradiance: dict,
        actual_weather_data: pd.Series) -> (
        eb_inputs.Inputs,
        eb_params.Params):
    raw_inputs = json_inputs.copy()

    raw_inputs.update(
        {"leaf_layers": vegetative_layers,
         "solar_inclination": actual_weather_data['solar_declination'],
         "wind_speed": actual_weather_data['wind_speed'],
         "vapor_pressure": actual_weather_data['vapor_pressure'],
         "vapor_pressure_deficit": actual_weather_data['vapor_pressure_deficit'],
         "air_temperature": actual_weather_data['air_temperature'],
         "incident_photosynthetically_active_radiation": {
             'direct': actual_weather_data['incident_direct_irradiance'],
             'diffuse': actual_weather_data['incident_diffuse_irradiance']},
         "absorbed_photosynthetically_active_radiation": absorbed_par_irradiance})

    return raw_inputs, json_params


def solve_energy_balance(
        correct_stability: bool,
        vegetative_layers: dict,
        leaf_class_type: str,
        absorbed_par_irradiance: dict,
        actual_weather_data: pd.Series) -> eb_solver.Solver:
    kwargs = {k: v for k, v in locals().items() if k not in ('leaf_class_type', 'correct_stability')}
    inputs_dict, params_dict = get_energy_balance_inputs_and_params(**kwargs)

    solver = eb_solver.Solver(leaves_category=leaf_class_type,
                              inputs_dict=inputs_dict,
                              params_dict=params_dict)
    solver.run(is_stability_considered=correct_stability)

    return solver


def get_variable(
        var_to_get: str,
        one_step_solver: eb_solver.Solver,
        leaf_class_type: str) -> dict:
    if leaf_class_type == 'lumped':
        res = {index: {'lumped': getattr(layer, var_to_get)} for index, layer in one_step_solver.crop.items()}
    else:
        res = {index: {'sunlit': getattr(layer['sunlit'], var_to_get), 'shaded': getattr(layer['shaded'], var_to_get)}
               for index, layer in one_step_solver.crop.items() if index != -1}
        res.update({-1: {'lumped': getattr(one_step_solver.crop[-1], var_to_get)}})
    return res


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
        ax.grid()
        ax.set(xlabel='hour', ylim=(-10, 50)),
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

    for component_index in component_indexes:
        abs_irradiance = summary_data[component_index]
        if component_index == -1:
            ax.plot(range(24), abs_irradiance, label=f'soil', color='k', linewidth=2)
        else:
            if leaf_class == 'lumped':
                ax.plot(range(24), summary_data[component_index], label=f'abs {component_index}')
            else:
                ax.plot(range(24), summary_data[component_index]['sunlit'], label=f'abs sunlit {component_index}')
                ax.plot(range(24), summary_data[component_index]['shaded'], label=f'abs shaded {component_index}')
    pass


def plot_temperature_dynamics(ax: plt.axis,
                              temperature_air: pd.Series,
                              simulation_case: str,
                              all_cases_data: dict):
    summary_data = get_summary_data(simulation_case, all_cases_data)
    _, leaf_class = simulation_case.split('_')
    component_indexes = summary_data.keys()

    ax.set_title(simulation_case.replace('_', ' '))

    hours = range(24)
    ax.plot(hours, temperature_air, label='air', color='k', linestyle='--', linewidth=2)

    for component_index in component_indexes:
        component_temperature = summary_data[component_index]
        if component_index == -1:
            ax.plot(hours, [(v - 273.15 if v > 273.15 else None) for v in component_temperature],
                    label=f'soil', color='k', linewidth=2)
        else:
            if leaf_class == 'lumped':
                ax.plot(hours, [(v - 273.15 if v > 273.15 else None) for v in summary_data[component_index]],
                        label=f'{component_index}')
            else:
                ax.plot(hours,
                        [(v - 273.15 if v > 273.15 else None) for v in summary_data[component_index]['sunlit']],
                        label=f'sunlit {component_index}')
                ax.plot(hours,
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


def plot_canopy_variable(
        all_cases_solver: eb_canopy,
        variable_to_plot: str,
        axes: plt.axis = None,
        return_axes: bool = False):
    hours = range(24)
    cases = all_cases_solver.keys()

    if axes is None:
        _, axes = plt.subplots(ncols=len(cases), sharex='all', sharey='all', figsize=(15, 5))

    for i, case in enumerate(cases):
        ax = axes[i]
        ax.plot(hours, [getattr(all_cases_solver[case][h].crop.state_variables, variable_to_plot) for h in hours])
    if return_axes:
        return axes
    else:
        [ax.set_xlabel('hours') for ax in axes[:]]
        plt.suptitle(variable_to_plot)
        plt.savefig(f'figs/{variable_to_plot}.png')
        plt.close()


def plot_energy_balance_components(
        h_solver: list,
        variable_to_plot: str,
        ax: plt.Subplot,
        return_ax: bool = False):
    hours = range(24)

    if variable_to_plot == 'soil_heat_flux':
        ax.plot(hours, [getattr(h_solver[h].crop[-1], 'heat_flux') for h in hours],
                label=UNITS_MAP[variable_to_plot][0])
    elif variable_to_plot == 'energy_balance':
        ax.plot(hours, [h_solver[h].energy_balance for h in hours], 'k--', label=UNITS_MAP[variable_to_plot][0])
    else:
        ax.plot(hours, [getattr(h_solver[h].crop.state_variables, variable_to_plot) for h in hours],
                label=UNITS_MAP[variable_to_plot][0])

    if return_ax:
        return ax
    else:
        ax.set_xlabel('hours')
        fig = ax.get_figure()
        fig.suptitle(variable_to_plot)
        fig.savefig(f'figs/{variable_to_plot}.png')
        plt.close()


def plot_energy_balance(solvers: dict):
    eb_components = [
        'net_radiation', 'sensible_heat_flux', 'total_penman_monteith_evaporative_energy', 'soil_heat_flux',
        'energy_balance']
    models = solvers.keys()
    fig, axes = plt.subplots(ncols=len(models), sharex='all', sharey='all', figsize=(15, 5))
    for model, ax in zip(models, axes):
        for eb_component in eb_components:
            ax = plot_energy_balance_components(h_solver=solvers[model], variable_to_plot=eb_component, ax=ax,
                                                return_ax=True)
        ax.grid()
        ax.legend()
    axes[0].set_ylabel(r'$\mathregular{Energy\/[W\/m^{-2}_{ground}]}$')
    fig.savefig(f'figs/energy_balance.png')
    pass


def plot_stability_terms(solvers: dict):
    # hourly plots
    terms = ['sensible_heat_flux', 'monin_obukhov_length', 'stability_correction_for_momentum',
             'stability_correction_for_heat']
    models = solvers.keys()
    fig, axes = plt.subplots(ncols=len(models), nrows=len(terms), sharex='all', sharey='row', figsize=(15, 5))

    for i, term in enumerate(terms):
        plot_canopy_variable(all_cases_solver=solvers, variable_to_plot=term, axes=axes[i, :], return_axes=True)
        axes[i, 0].set_ylabel(' '.join(UNITS_MAP[term]))

    [ax.set_xlabel('hour') for ax in axes[-1, :]]
    fig.savefig(f'figs/stability_terms.png')
    plt.close('all')


def plot_universal_functions(solvers, measurement_height: float = 2):
    # universal functions plots
    x = []
    phi_h = []
    phi_m = []
    for m_solvers in solvers.values():
        for h_solver in m_solvers:
            state = h_solver.crop.state_variables
            x.append((measurement_height - state.zero_displacement_height) / state.monin_obukhov_length)
            phi_h.append(state.stability_correction_for_heat)
            phi_m.append(state.stability_correction_for_momentum)
    _, axes = plt.subplots(ncols=2, sharex='all', sharey='all')
    axes[0].scatter(x, phi_h)
    axes[0].set_ylabel(' '.join(UNITS_MAP['stability_correction_for_heat']))
    axes[1].scatter(x, phi_m)
    axes[1].set_ylabel(' '.join(UNITS_MAP['stability_correction_for_momentum']))
    for ax in axes:
        ax.set_xlabel(r'$\mathregular{\frac{z_m-d}{L}\/[m]}$')
        ax.grid()
    plt.savefig(f'figs/universal_functions.png')
    plt.close('all')


if __name__ == '__main__':
    (Path(__file__).parent / 'figs').mkdir(exist_ok=True)
    weather_data = foo.get_weather_data()
    correct_for_stability = False
    irradiance = {}
    temperature = {}
    layers = {}
    solver_group = {}

    for canopy_type, leaves_type in (('bigleaf', 'lumped'),
                                     ('bigleaf', 'sunlit-shaded'),
                                     ('layered', 'lumped'),
                                     ('layered', 'sunlit-shaded')):
        print('-' * 50)
        print(f"{canopy_type} - {leaves_type}")
        canopy_layers = {0: sum(leaf_layers.values())} if canopy_type == 'bigleaf' else leaf_layers
        layers.update({f'{canopy_type} {leaves_type}': canopy_layers})
        hourly_absorbed_irradiance = []
        hourly_temperature = []
        hourly_solver = []

        for date, w_data in weather_data.iterrows():
            incident_direct_irradiance = w_data['incident_direct_irradiance']
            incident_diffuse_irradiance = w_data['incident_diffuse_irradiance']
            solar_inclination = w_data['solar_declination']
            wind_speed = w_data['wind_speed']
            vapor_pressure_deficit = w_data['vapor_pressure_deficit']
            vapor_pressure = w_data['vapor_pressure']
            air_temperature = w_data['air_temperature']
            print(date)
            absorbed_irradiance = calc_absorbed_irradiance(
                is_bigleaf=(canopy_type == 'bigleaf'),
                is_lumped=(leaves_type == 'lumped'),
                incident_direct_par_irradiance=incident_direct_irradiance,
                incident_diffuse_par_irradiance=incident_diffuse_irradiance,
                solar_inclination_angle=solar_inclination)

            hourly_absorbed_irradiance.append(absorbed_irradiance)

            energy_balance_solver = solve_energy_balance(
                vegetative_layers=canopy_layers,
                leaf_class_type=leaves_type,
                absorbed_par_irradiance=absorbed_irradiance,
                actual_weather_data=w_data,
                correct_stability=correct_for_stability)

            hourly_solver.append(energy_balance_solver)
            hourly_temperature.append(get_variable(
                var_to_get='temperature',
                one_step_solver=energy_balance_solver,
                leaf_class_type=leaves_type))

        irradiance[f'{canopy_type}_{leaves_type}'] = hourly_absorbed_irradiance
        temperature[f'{canopy_type}_{leaves_type}'] = hourly_temperature
        solver_group[f'{canopy_type}_{leaves_type}'] = hourly_solver

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

    plot_canopy_variable(
        all_cases_solver=solver_group,
        variable_to_plot='source_temperature')

    plot_energy_balance(
        solvers=solver_group)

    plot_stability_terms(solvers=solver_group)

    if correct_for_stability:
        plot_universal_functions(solvers=solver_group, measurement_height=2)
