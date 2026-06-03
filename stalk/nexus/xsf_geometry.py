#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path

from nexus import Structure

from stalk.nexus.nexus_structure import NexusStructure
from stalk.io.geometry_writer import GeometryWriter
from stalk.io.geometry_loader import GeometryLoader
from stalk.params.geometry_result import GeometryResult


class XsfGeometry(GeometryLoader, GeometryWriter):
    _suffix = 'structure.xsf'

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        scale=1.0,
        c_pos=None,
        **kwargs
    ):
        GeometryLoader.__init__(self, args=args, scale=scale, c_pos=c_pos, **kwargs)
        # GeometryWriter needs no additional initialization
    # end def

    def _load(self, filename) -> GeometryResult:
        # Using Nexus implementation to load XSF
        s = Structure()
        filename = self.get_filename(filename)
        if not filename.is_file():
            raise FileNotFoundError(f'File {filename} not found. Cannot load geometry.')
        # end if
        s.read_xsf(str(filename))
        return GeometryResult(s.pos, axes=s.axes, elem=s.elem)
    # end def

    def _write(
        self,
        structure: NexusStructure,
        path: Path,
        **kwargs,
    ):
        if not isinstance(structure, NexusStructure):
            raise TypeError('Presently only NexusStructure can be written to XSF file. Aborting.')
        # end ifs
        s = structure.get_nexus_structure()
        s.write_xsf(path)
    # end def

# end class
