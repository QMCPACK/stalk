#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path

from numpy import array

from stalk.params.parameter_set import ParameterSet
from stalk.params.parameter_structure import ParameterStructure
from stalk.io.geometry_writer import GeometryWriter
from stalk.io.geometry_loader import GeometryLoader
from stalk.params.geometry_result import GeometryResult


class XyzGeometry(GeometryWriter, GeometryLoader):
    _suffix = 'structure.xyz'

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        **kwargs  # scale=1.0, c_pos=None, etc.
    ):
        GeometryLoader.__init__(self, args=args, **kwargs)
        # GeometryWriter needs no additional initialization
    # end def

    # Override loading hook
    def _load(self, filename) -> GeometryResult:
        data = self.load_result(filename, dtype=str, unpack=True, skiprows=2, rescale=False)
        el, x, y, z = data
        # Apply scaling only after type conversion
        pos = array([x, y, z], dtype=float).T / self.scale
        return GeometryResult(pos, axes=None, elem=el)
    # end def

    # Override writing hook
    def _write(self, structure: ParameterSet, path: Path, **kwargs):
        output = []
        if isinstance(structure, ParameterStructure):
            pos = structure.pos.copy()
            elem = structure.elem
        elif isinstance(structure, ParameterSet):
            pos = structure.params.copy()
            elem = 'p'
        else:
            raise TypeError(f'Cannot write to XYZ file: {structure}')
        # end if

        header = str(len(elem)) + '\n'
        fmt = '{:< 10f}'
        for el, pr in zip(elem, pos):
            row = [el]
            for p in pr:
                row.append(fmt.format(p))
            # end for
            output.append(row)
        # end for
        # Using TxtResult methods
        self.save_result(
            path,
            data=array(output),
            header=header,
            fmt='%s',
            comments=''
        )
    # end def

# end class
