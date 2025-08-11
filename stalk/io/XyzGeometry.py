#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import loadtxt, savetxt, array

from stalk.params.ParameterSet import ParameterSet
from stalk.params.ParameterStructure import ParameterStructure
from stalk.io.GeometryWriter import GeometryWriter
from stalk.io.GeometryLoader import GeometryLoader
from stalk.params.GeometryResult import GeometryResult
from stalk.util.util import PL


class XyzGeometry(GeometryLoader, GeometryWriter):

    def __init__(
        self,
        args={}
    ):
        self.args = args
    # end def

    def _load(self, path, suffix='structure.xyz', c_pos=1.0):
        el, x, y, z = loadtxt(
            PL.format(path, suffix),
            dtype=str,
            unpack=True,
            skiprows=2
        )
        pos = array([x, y, z], dtype=float).T * c_pos
        return GeometryResult(pos, axes=None, elem=el)
    # end def

    def __write__(self, structure, path, suffix='structure.xyz', c_pos=1.0):
        output = []
        if isinstance(structure, ParameterStructure):
            pos = structure.pos
            elem = structure.elem
        elif isinstance(structure, ParameterSet):
            pos = structure.params
            elem = 'p'
        else:
            raise TypeError(f'Cannot write to XYZ file: {structure}')
        # end if

        header = str(len(elem)) + '\n'
        fmt = '{:< 10f}'
        for el, pr in zip(elem, pos * c_pos):
            row = [el]
            for p in pr:
                row.append(fmt.format(p))
            # end for
            output.append(row)
        # end for
        savetxt(
            PL.format(path, suffix),
            array(output),
            header=header,
            fmt='%s',
            comments=''
        )
    # end def

# end class
