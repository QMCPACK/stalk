#!/usr/bin/env python3
"""Base class for representing a set of parameters to optimize"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array, isscalar, random, ndarray, linspace
from copy import deepcopy

from stalk.params.linesearch_point import LineSearchPoint
from stalk.params.parameter import Parameter


class ParameterSet(LineSearchPoint):
    _param_list: list[Parameter] = []
    _samples = None  # samples for effective variance estimation
    label = ''  # label for identification
    file_path = None  # field to be used in file I/O mode

    def __init__(
        self,
        params=None,  # List of scalars or Parameter objects
        params_err=None,
        value=None,
        error=0.0,
        label=None,
    ):
        self.label = label
        self.params_list = params
        self.params_err = params_err
        if value is not None:
            self.value = value
            self.error = error
        # end if
    # end def

    @property
    def params_list(self):
        return [p for p in self._param_list if isinstance(p, Parameter)]
    # end def

    @params_list.setter
    def params_list(self, params_list=[]):
        if params_list is None:
            params_list = []
        # end if
        p_list = []
        for p, param in enumerate(params_list):
            if isinstance(param, Parameter):
                parameter = param
            elif isscalar(param):
                label = f'p{p}'
                parameter = Parameter(param, 0.0, unit=None, label=label)
            else:
                raise ValueError('Parameter is unsupported type: ' + str(param))
            # end if
            p_list.append(parameter)
        # end for
        self._param_list = p_list
        # Reset value upon params init
        self.reset_value()
    # end def

    @property
    def params(self):
        if len(self) > 0:
            return array([p.value for p in self.params_list])
        # end if
    # end def

    @params.setter
    def params(self, params: ndarray | list[Parameter] | list):
        if len(self) == 0:
            # Try to initialize parameters with default values
            self.params_list = params
        elif len(params) != len(self):
            raise ValueError(f'Inconsistent size of new params {len(params)} vs {len(self)}')
        else:
            for p, param in zip(params, self.params_list):
                param.value = p
                # Reset parameter error upon parameter change
                param.error = 0.0
            # end if
            self.reset_value()
        # end if
    # end def

    @property
    def params_err(self):
        if len(self) > 0:
            return array([p.error for p in self.params_list])
        # end if
    # end def

    @params_err.setter
    def params_err(self, params_err: ndarray | list):
        if params_err is None:
            params_err = len(self) * [0.0]
        elif len(params_err) != len(self):
            raise ValueError(f'Inconsistent size of new params {len(params_err)} vs {len(self)}')
        # end if
        for p, param in zip(params_err, self.params_list):
            param.error = p
        # end if
    # end def

    @property
    def samples(self):
        return self._samples
    # end def

    @samples.setter
    def samples(self, samples):
        if samples is None:
            self._samples = None
        elif int(samples) > 0:
            self._samples = int(samples)
        else:
            raise ValueError(f'Samples must be > 0 integer, provided: {samples}')
        # end if
    # end def

    def shift_params(self, shifts):
        if len(shifts) != len(self):
            raise ValueError('Shifts has wrong dimensions!')
        # end if
        for param, shift in zip(self.params_list, shifts):
            param.shift(shift)
        # end for
        self.reset_value()
    # end def

    def copy(
        self,
        params=None,
        params_err=None,
        label=None,
        offset=None
    ):
        paramset = deepcopy(self)
        if offset is not None:
            paramset.offset = offset
        # end if
        if params is not None:
            paramset.params = params
        # end if
        if params_err is not None:
            paramset.params_err = params_err
        # end if
        if label is not None:
            paramset.label = label
        # end if
        return paramset
    # end def

    def get_params_distribution(self, N=100):
        return [self.params + self.params_err * g for g in random.randn(N, len(self))]
    # end def

    def check_consistency(self):
        return True
    # end def

    def __sub__(self, other):
        if isinstance(other, ParameterSet) and len(self) > 0 and len(other) == len(self):
            result = self.copy()
            result.shift_params(-other.params)
            return result
        else:
            raise ValueError(f'Cannot subtract {repr(other)} from {repr(self)}')
        # end if
    # end def

    # Mean squared distance between this and other ParameterSet
    def distance2(self, other):
        diff = self - other
        return sum(array(diff.params)**2)
    # end def

    # Mean unsigned distance between this and other ParameterSet
    def distance(self, other):
        return self.distance2(other)**0.5
    # end def

    def __str__(self):
        string = self.__class__.__name__
        if self.label is not None or self.label != '':
            string += ' ({})'.format(self.label)
        # end if
        if self.params is None:
            string += '\n  params: not set'
        else:
            string += '\n  params:'
            for param in self.params_list:
                string += '\n    ' + str(param)
            # end for
        # end if
        if self.value is not None:
            string += f'\n  value: {self.value}'
            if self.error > 0:
                string += f' +/- {self.error}'
            # end if
        # end if
        return string
    # end def

    def __len__(self):
        return len(self.params_list)
    # end def

# end class


def interpolate_params(structure_a: ParameterSet, structure_b: ParameterSet, num_int):
    scales = linspace(0.0, 1.0, num_int + 2)
    dparams = structure_b.params - structure_a.params
    traj = []
    for scale in scales:
        new_params = structure_a.params + scale * dparams
        structure = structure_a.copy(params=new_params)
        traj.append(structure)
    # end for
    return traj
# end def
