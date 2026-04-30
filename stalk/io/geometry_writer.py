#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from pathlib import Path

from stalk.util import directorize
from stalk.util.args_container import ArgsContainer


class GeometryWriter(ArgsContainer):

    def write(self, structure, path='', **kwargs):
        '''The Geometry writer must accept a "structure" and a "path" to output file
        '''
        # Hot update of args
        args = self.get_updated(kwargs)

        # If suffix is present, write to 'path/suffix', else write to 'path'
        suffix = args.pop('suffix', None)
        if suffix is None:
            filename = path
        else:
            filename = f'{directorize(path)}{suffix}'
        # end if
        if Path(filename).is_dir():
            raise IsADirectoryError(f"The target file {filename} is a directory!")
        # end if
        self._write(structure, filename, **args)
    # end def

    # The actual writing function must be overridden
    def _write(self, structure, filename, **kwargs):
        raise NotImplementedError("Implement _write(structure, filename) function in inherited class.")
    # end def

# end class
