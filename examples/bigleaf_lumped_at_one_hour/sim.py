from pathlib import Path

from matplotlib import pyplot

from crop_energy_balance.crop import Crop
from crop_energy_balance.formalisms.weather import convert_kelvin_to_celsius
from crop_energy_balance.solver import Solver


def plot_temperature_profile(canopy_object: Crop, fig_path: Path):
    soil_component_key = -1
    leaf_layer_key = 0
    soil_component_temperature = convert_kelvin_to_celsius(canopy_object[-1].temperature)
    leaf_component_temperature = convert_kelvin_to_celsius(canopy_object[0].temperature)

    fig, ax = pyplot.subplots()
    ax.grid(zorder=0)
    ax.scatter(soil_component_temperature, soil_component_key,
               marker='o', c='brown', edgecolors=None, label='soil', zorder=3)
    ax.scatter(leaf_component_temperature, leaf_layer_key,
               marker='o', c='orange', edgecolors=None, label='lumped leaves', zorder=3)
    ax.axvline(convert_kelvin_to_celsius(canopy_object.inputs.air_temperature),
               label='air', zorder=3)
    ax.set(xlabel='temperature [°C]', ylabel='Component index [-]', xlim=(20, 30))
    ax.yaxis.set_ticks([-1, 0])
    ax.legend()
    fig.savefig(fig_path)
    pyplot.close()


if __name__ == '__main__':
    root_pth = Path(__file__).parent

    solver = Solver(leaves_category='lumped',
                    inputs_path=root_pth / 'inputs_well_watered.json',
                    params_path=root_pth / 'params.json')
    solver.run()

    plot_temperature_profile(solver.crop, root_pth / 'temperature_profile.png')
    print(f'*** iterations = {solver.iterations_number} **')
