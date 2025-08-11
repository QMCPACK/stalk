#!/usr/bin/env python3
"""Base class for representing a set of parameters to optimize"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array, isscalar, random
from copy import deepcopy

from stalk.params.LineSearchPoint import LineSearchPoint
from stalk.params.Parameter import Parameter


class ParameterSet(LineSearchPoint):
    _param_list: list[Parameter] = []
    label = ''  # label for identification
    file_path = None  # field to be used in file I/O mode

    @property
    def params_list(self):
        return [p for p in self._param_list if isinstance(p, Parameter)]
    # end def

    @property
    def params(self):
        if len(self) > 0:
            return array([p.value for p in self.params_list])
        # end if
    # end def

    @property
    def params_err(self):
        if len(self) > 0:
            return array([p.error for p in self.params_list])
        # end if
    # end def

    def __init__(
        self,
        params=None,  # List of scalars or Parameter objects
        params_err=None,
        units=None,
        value=None,
        error=0.0,
        label=None,
        labels=None
    ):
        self.label = label
        if params is not None:
            self.init_params(
                params,
                errors=params_err,
                units=units,
                labels=labels
            )
        # end if
        if value is not None:
            self.value = value
            self.error = error
        # end if
    # end def

    def init_params(self, params, errors=None, units=None, labels=None):
        if errors is None:
            errors = len(params) * [errors]
        else:
            assert len(errors) == len(params)
        # end if
        if units is None or isinstance(units, str):
            units = len(params) * [units]
        else:
            assert len(units) == len(params)
        # end if
        if labels is None:
            labels = len(params) * [labels]
        else:
            assert len(labels) == len(labels)
        # end if
        p_list = []
        for p, (param, error, unit, label) in enumerate(zip(params, errors, units, labels)):
            if isinstance(param, Parameter):
                parameter = param
            elif isscalar(param):
                lab = label if label is not None else 'p{}'.format(p)
                parameter = Parameter(param, error, unit=unit, label=lab)
            else:
                raise ValueError('Parameter is unsupported type: ' + str(param))
            # end if
            p_list.append(parameter)
        # end for
        self._param_list = p_list
    # end def

    def set_params(self, params, params_err=None):
        # If params have not been initiated yet, do it now without extra info
        if len(self) > 0:
            assert len(params) == len(self)
        else:
            self.init_params(params, params_err)
        # end if
        if params_err is None:
            params_err = len(params) * [0.0]
        # end if
        for sparam, param, error in zip(self.params_list, params, params_err):
            assert isinstance(sparam, Parameter)
            sparam.value = param
            sparam.error = error
        # end for
        self.reset_value()
    # end def

    def shift_params(self, shifts):
        if len(shifts) != len(self):
            raise ValueError('Shifts has wrong dimensions!')
        # end if
        for param, shift in zip(self.params_list, shifts):
            assert isinstance(param, Parameter)
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
            paramset.set_params(params, params_err)
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
        if isinstance(other, ParameterSet) and len(other) == len(self):
            result = self.copy()
            result.shift_params(-other.params)
        else:
            raise ValueError(f'Cannot subtract {repr(other)} from {repr(self)}')
        # end if
    # end def

    def __str__(self):
        string = self.__class__.__name__
        if self.label is not None:
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
        return string
    # end def

    def __len__(self):
        return len(self.params_list)
    # end def

# end class
