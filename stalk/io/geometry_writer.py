#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from stalk.io.txt_data import TxtData
from stalk.util.args_container import ArgsContainer


class GeometryWriter(ArgsContainer, TxtData):

    def write(self, structure, path: str):
        # Writing hook
        self._write(structure, path, **self.args)
    # end def

    # The actual writing function must be overridden in a derived class
    def _write(self, structure, filename, **kwargs):
        raise NotImplementedError("Implement _write(structure, filename) function in inherited class.")
    # end def

# end class
