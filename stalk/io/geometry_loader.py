#!/usr/bin/env python3

import warnings
from numpy import nan

from numpy import isscalar

from stalk.params.geometry_result import GeometryResult
from stalk.params.parameter_set import ParameterSet
from stalk.util.args_container import ArgsContainer
from stalk.util.util import check_result_file

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class GeometryLoader(ArgsContainer):

    def load(self, path, c_pos=None, only_warn=False, **kwargs) -> GeometryResult:
        '''The Geometry loader must accept a "path" to input file and return GeometryResult.
        '''
        # Hot update of args
        args = self.get_updated(kwargs)
        # Keep for backward compatibility
        if isscalar(c_pos):
            args['scale'] = c_pos**-1
        # end if
        scale = args.pop('scale', 1.0)

        try:
            filename = check_result_file(path, args)
        except FileNotFoundError as e:
            if only_warn:
                warnings.warn(repr(e))
                return GeometryResult(nan)
            else:
                raise e
            # end if
        # end try
        res = self._load(filename, **args)
        # Rescale to model units
        res.rescale(scale)
        print(f'Loaded geometry from {filename}')
        return res
    # end def

    def load_or_relax(
        self,
        path,
        relax_func: callable,
        structure: ParameterSet,
        **kwargs  # relax kwargs
    ) -> ParameterSet:
        try:
            res = self.load(path)
            if not callable(relax_func):
                raise TypeError('The relax_func must be callable and return GeometryResult!')
            # end if
        except FileNotFoundError:
            # Try to relax
            relax_func(structure.copy(), **kwargs)
            # Then, try to load again
            res = self.load(path)
        # end try
        return structure.copy(pos=res.pos, axes=res.axes)
    # end def

    # The actual loading function must be overridden and return a GeometryResult object
    def _load(self, filename, **kwargs) -> GeometryResult:
        raise NotImplementedError("Implement _load(filename) function in inherited class.")
    # end def

# end class
