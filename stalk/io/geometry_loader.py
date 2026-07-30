#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.io.txt_data import TxtData
from stalk.params.geometry_result import GeometryResult
from stalk.params.parameter_set import ParameterSet
from stalk.util.args_container import ArgsContainer


class GeometryLoader(ArgsContainer, TxtData):
    _suffix = 'structure.dat'

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        scale=1.0,
        c_pos=None,
        **kwargs,
    ):
        # Backward compatibility: if c_pos is given, it overrides scale
        if c_pos is not None:
            scale = c_pos**-1
        # end if
        args.update(**kwargs)
        # suffix = None means we'll use class-level default
        suffix = args.pop('suffix', None)
        TxtData.__init__(self, suffix=suffix, scale=scale)
        ArgsContainer.__init__(self, **args)
    # end def

    def load(self, path) -> GeometryResult:
        # Loading hook
        res = self._load(path, **self.args)
        print(f'Loaded geometry from {self.get_filename(path)}.')
        return res
    # end def

    def load_or_relax(
        self,
        path,
        relax_func: callable,
        structure: ParameterSet,
        label='relax',
        **kwargs  # relax kwargs
    ) -> ParameterSet:
        try:
            res = self.load(path)
            if not callable(relax_func):
                raise TypeError('The relax_func must be callable and write the geometry result file to the same path')
            # end if
        except FileNotFoundError:
            # Try to relax
            relax_func(structure.copy(), **kwargs)
            # Then, try to load again
            res = self.load(path)
        # end try
        return structure.copy(pos=res.pos, axes=res.axes, label=label)
    # end def

    # The actual loading function must be overridden and return a GeometryResult object
    def _load(self, path: str, **kwargs) -> GeometryResult:
        raise NotImplementedError("Implement _load(filename) function in inherited class.")
    # end def

# end class
