from pathlib import Path
from matplotlib import pyplot

from crop_energy_balance.crop import Canopy
from crop_energy_balance.inputs import SunlitShadedInputs
from crop_energy_balance.params import SunlitShadedParams
from crop_energy_balance.solver import Solver
from crop_energy_balance.utils import convert_kelvin_to_celsius


def plot_temperature_profile(canopy_object: Canopy, fig_path: Path):
    soil_component_key = -1
    leaf_layer_key = 0
    soil_component_temperature = convert_kelvin_to_celsius(canopy_object[-1].temperature)
    sunlit_component_temperature = convert_kelvin_to_celsius(canopy_object[leaf_layer_key]['sunlit'].temperature)
    shaded_component_temperature = convert_kelvin_to_celsius(canopy_object[leaf_layer_key]['shaded'].temperature)

    fig, ax = pyplot.subplots()
    ax.grid(True, zorder=0)
    ax.scatter(soil_component_temperature, soil_component_key,
               marker='o', c='brown', edgecolors=None, label='soil', zorder=3)
    ax.scatter(sunlit_component_temperature, leaf_layer_key,
               marker='o', c='yellow', edgecolors=None, label='sunlit leaves', zorder=3)
    ax.scatter(shaded_component_temperature, leaf_layer_key,
               marker='o', c='orange', edgecolors=None, label='shaded leaves', zorder=3)
    ax.axvline(convert_kelvin_to_celsius(canopy_object.inputs.air_temperature), label='air')
    ax.set(xlabel='temperature [°C]', ylabel='Component index [-]', xlim=(20, 30))
    ax.yaxis.set_ticks([-1, 0])
    ax.legend()
    fig.savefig(fig_path)
    pyplot.close()


if __name__ == '__main__':
    root_pth = Path(__file__).parent

    inputs = SunlitShadedInputs(inputs_path=root_pth / 'inputs_well_watered.json')
    params = SunlitShadedParams(params_path=root_pth / 'params.json')
    params.update(inputs=inputs)

    canopy = Canopy(leaves_category='sunlit-shaded', inputs=inputs, params=params)
    solver = Solver(canopy=canopy, inputs=inputs, params=params)
    solver.run()

    plot_temperature_profile(canopy, root_pth / 'temperature_profile.png')
    print(f'*** iterations = {solver.iterations_number} **')
