#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from nexus import Structure

from stalk.nexus.nexus_structure import NexusStructure
from stalk.io.geometry_writer import GeometryWriter
from stalk.io.geometry_loader import GeometryLoader
from stalk.params.geometry_result import GeometryResult


class XsfGeometry(GeometryLoader, GeometryWriter):

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        suffix='structure.xsf',
        **kwargs
    ):
        my_args = {'suffix': suffix}
        my_args.update(**args, **kwargs)
        super().__init__(**my_args)
    # end def

    def _load(self, filename) -> GeometryResult:
        # Using Nexus implementation to load XSF
        s = Structure()
        s.read_xsf(filename)
        return GeometryResult(s.pos, axes=s.axes, elem=s.elem)
    # end def

    def _write(
        self,
        structure: NexusStructure,
        filename: str,
        **kwargs,
    ):
        if not isinstance(structure, NexusStructure):
            raise TypeError('Presently only NexusStructure can be written to XSF file. Aborting.')
        # end ifs
        s = structure.get_nexus_structure()
        s.write_xsf(filename)
    # end def

# end class
