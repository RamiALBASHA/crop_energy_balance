from pathlib import Path

from crop_energy_balance.crop import Canopy
from crop_energy_balance.inputs import Inputs
from crop_energy_balance.params import Params
from crop_energy_balance.solver import Solver
from crop_energy_balance.formalisms.weather import convert_kelvin_to_celsius
from matplotlib import pyplot


def plot_temperature_profile(canopy_object: Canopy, fig_path: Path):
    soil_component_key = -1
    leaf_layer_keys = list(canopy_object.inputs.leaf_layers.keys())
    soil_component_temperature = convert_kelvin_to_celsius(canopy_object[-1].temperature)
    sunlit_component_temperature = [convert_kelvin_to_celsius(canopy_object[k]['sunlit'].temperature)
                                    for k in leaf_layer_keys]
    shaded_component_temperature = [convert_kelvin_to_celsius(canopy_object[k]['shaded'].temperature)
                                    for k in leaf_layer_keys]

    fig, ax = pyplot.subplots()
    ax.grid(zorder=0)
    ax.scatter(soil_component_temperature, soil_component_key,
               marker='o', c='brown', edgecolors=None, label='soil', zorder=3)
    ax.plot(sunlit_component_temperature, leaf_layer_keys,
            marker=None, c='yellow', label='sunlit leaves', zorder=3)
    ax.plot(shaded_component_temperature, leaf_layer_keys,
            marker=None, c='orange', label='shaded leaves', zorder=3)
    ax.axvline(convert_kelvin_to_celsius(canopy_object.inputs.air_temperature), label='air')
    ax.set(xlabel='temperature [°C]', ylabel='Component index [-]', xlim=(20, 30))
    ax.legend()
    fig.savefig(fig_path)
    pyplot.close()


if __name__ == '__main__':
    root_pth = Path(__file__).parent

    inputs = Inputs(inputs_path=root_pth / 'inputs_well_watered.json')
    params = Params(params_path=root_pth / 'params.json')
    params.update(inputs=inputs)

    canopy = Canopy(leaves_category='sunlit-shaded', inputs=inputs, params=params)
    solver = Solver(canopy=canopy, inputs=inputs, params=params)
    solver.run()

    plot_temperature_profile(canopy, root_pth / 'temperature_profile.png')
    print(f'*** iterations = {solver.iterations_number} **')
