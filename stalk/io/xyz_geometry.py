#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import loadtxt, savetxt, array

from stalk.params.parameter_set import ParameterSet
from stalk.params.parameter_structure import ParameterStructure
from stalk.io.geometry_writer import GeometryWriter
from stalk.io.geometry_loader import GeometryLoader
from stalk.params.geometry_result import GeometryResult


class XyzGeometry(GeometryLoader, GeometryWriter):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='structure.xyz',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(self, filename) -> GeometryResult:
        el, x, y, z = loadtxt(
            filename,
            dtype=str,
            unpack=True,
            skiprows=2
        )
        pos = array([x, y, z], dtype=float).T
        return GeometryResult(pos, axes=None, elem=el)
    # end def

    def _write(self, structure: ParameterSet, filename):
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
        savetxt(
            filename,
            array(output),
            header=header,
            fmt='%s',
            comments=''
        )
    # end def

# end class
