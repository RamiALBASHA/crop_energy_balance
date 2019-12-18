from crop_energy_balance.inputs import LumpedInputs, SunLitShadedInputs
from crop_energy_balance.params import Params


class LeafLayer:
    def __init__(self,
                 index: int,
                 upper_cumulative_leaf_area_index: float,
                 thickness: float):
        self.index = index
        self.upper_cumulative_leaf_area_index = upper_cumulative_leaf_area_index
        self.thickness = thickness
        self.absorbed_irradiance = {}

    def calc_absorbed_irradiance(self,
                                 inputs: LumpedInputs or SunlitShadedInputs,
                                 params: Params):
        pass


class SunlitShadedLeafLayer(LeafLayer):
    def __init__(self, **kwargs):
        LeafLayer.__init__(self, **kwargs)


class Shoot(dict):
    def __init__(self,
                 leaves_category: str,
                 inputs: LumpedInputs or SunlitShadedInputs,
                 params: Params):
        """Creates a class:`Shoot` object having either 'lumped' leaves or 'sunlit-shaded' leaves.
        Args:
            leaves_category: one of ('lumped', 'sunlit-shaded')
            inputs: see class`LumpedInputs` and `SunlitShadedInputs`
            params: see class`Params`
        Notes:
            The created shoot can implicitly be 'big-leaf' or a 'layered'. If the attribute `leaf_layers` of the
                :Class:`inputs` object has only one layer, then the resulting shoot is a 'big-leaf', otherwise if the
                dictionary has more than one key, then the shoot is 'layered'.
            Leaf layers indexes in `leaf_layers` must be ordered so that the youngest leaf layer has the highest index
                value, and inversely, the oldest leaf layer has the least value. Not respecting this order will
                definitely lead to erroneous calculations.
        """

        dict.__init__(self)

        self.inputs = inputs
        self.params = params
        self._leaf_layer_indexes = list(reversed(sorted(inputs.leaf_layers.keys())))

        self.set_leaf_layers(leaves_category)
        self.calc_absorbed_irradiance()

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def set_leaf_layers(self, leaves_category: str):
        """Sets leaf layers of the shoot.
        Args:
            leaves_category: one of ('lumped', 'sunlit-shaded')
        """

        upper_cumulative_leaf_area_index = 0.0
        for index in self._leaf_layer_indexes:
            layer_thickness = self.inputs.leaf_layers[index]
            if leaves_category == 'lumped':
                self[index] = LumpedLeafLayer(index,
                                              upper_cumulative_leaf_area_index,
                                              layer_thickness)
            else:
                self[index] = SunlitShadedLeafLayer(index,
                                                    upper_cumulative_leaf_area_index,
                                                    layer_thickness,
                                                    self.params)

            upper_cumulative_leaf_area_index += layer_thickness

    def calc_absorbed_irradiance(self):
        """Calculates the absorbed irradiance by shoot's layers.
        """
        for index in self._leaf_layer_indexes:
            self[index].calc_absorbed_irradiance(self.inputs, self.params)
