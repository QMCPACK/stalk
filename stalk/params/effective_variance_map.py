#!/usr/bin/env python3
'''Class to map parameter subspaces to EffectiveVariance estimates.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.params.parameter_set import ParameterSet
from stalk.params.effective_variance import EffectiveVariance


# Map relative scaling of the effective variance:
# {params} -> [var_eff_sample] -> mean(var_eff)
class EffectiveVarianceMap():
    _params = None
    scaling_map: list[(ParameterSet, EffectiveVariance)] = None

    @property
    def params(self):
        return self._params
    # end def

    @params.setter
    def params(self, params):
        if not isinstance(params, ParameterSet):
            raise TypeError('The provided params must be a ParameterSet')
        # end if
        if self.params is None or len(self.params) == len(params):
            self._params = params
        else:
            raise ValueError('The provided ParameterSet does not match the initial ParameterSet. Aborting.')
        # end if
    # end def

    def __init__(self, params: ParameterSet, var_eff: EffectiveVariance = None):
        self.params = params
        self.scaling_map = []
        if var_eff is not None:
            self.add_var_eff(params, var_eff)
        # end if
    # end def

    def add_var_eff(self, params: ParameterSet, var_eff: EffectiveVariance, thr=1e-6):
        # Check consistency between paramsets
        if not len(params) == len(self.params):
            raise AssertionError('The provided ParameterSet does not match the initial ParameterSet. Aborting.')
        # end if
        # implement a naive search loop for existing param data (WIP)
        found = False
        for p_this, var_this in self.scaling_map:
            # Params match to mean-squared deviation
            if sum((params.params - p_this.params)**2) < thr:
                var_this += var_eff
                found = True
                break
            # end if
        # end for
        if not found:
            # Add another entry
            self.scaling_map.append((params, var_eff))
        # end if
    # end def

    def get_samples(self, params: ParameterSet, error, nodes=1):
        if len(self) > 0:
            min_diff = 1e8
            i_nearest = 0
            # Implement a simple nearest-neightbor search
            for i in range(len(self)):
                p_diff = params.distance(self.scaling_map[i][0])
                if p_diff < min_diff:
                    i_nearest = i
                    min_diff = p_diff
                # end if
            # end for
            var_eff: EffectiveVariance = self.scaling_map[i_nearest][1]
            samples = var_eff.get_samples(error=error, nodes=nodes)
            return samples
        else:
            raise AssertionError('EffectiveVarianceMap is empty. Add data to allow estimation of samples')
        # end if
    # end def

    def __len__(self):
        return len(self.scaling_map)
    # end def

# end class
