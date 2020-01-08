from pathlib import Path
from matplotlib import pyplot

from crop_energy_balance.crop import Canopy
from crop_energy_balance.inputs import LumpedInputs
from crop_energy_balance.params import Params
from crop_energy_balance.solver import Solver
from crop_energy_balance.utils import convert_kelvin_to_celsius


def plot_temperature_profile(canopy_object: Canopy, fig_path: Path):
    component_index, component_temperature = zip(*[(k, convert_kelvin_to_celsius(v.temperature))
                                                   for k, v in canopy_object.items()])
    fig, ax = pyplot.subplots()
    ax.plot(component_temperature, component_index, 'g-', label='canopy')
    ax.axvline(convert_kelvin_to_celsius(canopy_object.inputs.air_temperature), label='air')
    ax.set(xlabel='temperature [°C]', ylabel='Component index [-]', xlim=(20, 30))
    ax.grid(True)
    fig.savefig(fig_path)
    pyplot.close()


if __name__ == '__main__':
    root_pth = Path(__file__).parent

    inputs = LumpedInputs(root_pth / 'inputs_well_watered.json')
    params = Params(root_pth / 'params.json')
    params.update(inputs=inputs)

    canopy = Canopy(leaves_category='lumped',
                    inputs=inputs,
                    params=params)
    Solver(canopy=canopy,
           inputs=inputs,
           params=params).run()

    plot_temperature_profile(canopy, root_pth / 'temperature_profile.png')
