#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import isscalar

from stalk.params.ParameterSet import ParameterSet
from stalk.params.PesResult import PesResult
from stalk.util.FunctionCaller import FunctionCaller


class PesLoader(FunctionCaller):
    _scale = 1.0

    def __init__(
        self,
        func,
        args={},
        scale=1.0,
    ):
        super().__init__(func, args)
        self.scale = scale
    # end def

    @property
    def scale(self):
        return self._scale
    # end def

    @scale.setter
    def scale(self, scale):
        if isscalar(scale) and scale > 0:
            self._scale = scale
        else:
            raise TypeError(f"The scale must scalar and >0, provided: {scale}")
        # end if
    # end

    def load(self, arg: ParameterSet | str, sigma=0.0, **kwargs):
        if isinstance(arg, str):
            structure = ParameterSet()
            structure.file_path = arg
        else:
            structure = arg
        # end if
        args = self.args.copy()
        args.update(kwargs)
        res = self._load(structure, **args)
        # Rescale to model units
        res.rescale(self.scale)
        # If a non-zero, artificial errorbar is requested, add it to result
        res.add_sigma(sigma)
        return res
    # end def

    def _load(self, structure, **kwargs):
        res = self.func(structure, **kwargs)
        if not isinstance(res, PesResult):
            raise AssertionError('The _load method must return a PesResult.')
        # end if
        return res
    # end def

# end class
