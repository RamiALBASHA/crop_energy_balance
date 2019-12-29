from pathlib import Path
from crop_energy_balance.inputs import Inputs
from crop_energy_balance.params import Params
from crop_energy_balance.crop import Canopy
from crop_energy_balance.solver import Solver

if __name__ == '__main__':
    root_pth = Path(__file__).parent

    inputs = Inputs(root_pth / 'inputs_well_watered.json')
    params = Params(root_pth / 'params.json')
    params.update(inputs=inputs)

    canopy = Canopy(leaves_category='lumped',
                    inputs=inputs,
                    params=params)
    Solver(canopy=canopy,
           inputs=inputs,
           params=params).run()

    x = 0
