#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path


class StalkPath():
    """
    A class to handle the path for object serialization and loading.
    """
    _path: Path = None  # field to be used in file I/O mode; If None, no I/O

    def __init__(
        self,
        path: str | Path | None
    ):
        self.path = path
    # end def

    @property
    def path(self) -> Path | None:
        return self._path
    # end def

    @path.setter
    def path(self, path):
        if path is None:
            self._path = None
        elif isinstance(path, (str, Path)):
            self._path = Path(path)
        else:
            raise ValueError(f'Path must be str or Path or None, provided: {path}')
        # end if
    # end def

    # Override to allow conditional path appending
    def __truediv__(
        self,
        other: str | Path | None
    ):
        if self.path is None or other is None:
            return None
        else:
            return self.path / Path(other)
        # end if
    # end def

# end class
