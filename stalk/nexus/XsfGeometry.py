#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
from nexus import Structure

from stalk.nexus.NexusStructure import NexusStructure
from stalk.io.GeometryWriter import GeometryWriter
from stalk.io.GeometryLoader import GeometryLoader
from stalk.params.GeometryResult import GeometryResult
from stalk.util.util import PL


class XsfGeometry(GeometryLoader, GeometryWriter):

    def __init__(
        self,
        args={}
    ):
        self.args = args
    # end def

    def _load(
        self,
        path: str,
        suffix='structure.xsf',
        c_pos=1.0
    ):
        s = Structure()
        fpath = PL.format(path, suffix)
        if not Path(fpath).exists():
            raise FileNotFoundError(f'Structure file {fpath} does not exist.')
        # end if
        s.read_xsf(fpath)
        pos = s.pos * c_pos
        if s.axes is not None:
            axes = s.axes * c_pos
        # end if
        return GeometryResult(pos, axes=axes, elem=s.elem)
    # end def

    def __write__(
        self,
        structure: NexusStructure,
        path: str,
        suffix='structure.xsf',
        c_pos=1.0
    ):
        if not isinstance(structure, NexusStructure):
            raise TypeError('Presently only NexusStructure can be written to XSF file. Aborting.')
        # end ifs
        s = structure.get_nexus_structure()
        s.pos *= c_pos
        if s.axes is not None:
            s.axes *= c_pos
        # end if
        s.write_xsf(PL.format(path, suffix))
    # end def

# end class
