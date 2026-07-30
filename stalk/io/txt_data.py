#!/usr/bin env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
from numpy import ndarray, loadtxt, array, isscalar, savetxt, nan


# TxtResult is a base class that contains suffix information and loading functionality for
# text-based result files, such as energy.dat
class TxtData:
    # Suffix is the filename in a path. The default suffix is preset in each derived class.
    _suffix: str = 'data.dat'
    _scale: float = None

    def __init__(self, suffix=None, scale=1.0):
        if suffix is not None:
            self.suffix = suffix
        # end if
        self.scale = scale
    # end def

    @property
    def scale(self):
        return self._scale
    # end def

    @scale.setter
    def scale(self, value):
        if isscalar(value) and value != 0.0:
            self._scale = value
        else:
            raise ValueError('Scale must be a non-zero scalar value.')
        # end if
    # end def

    @property
    def suffix(self):
        return self._suffix
    # end def

    @suffix.setter
    def suffix(self, value):
        if isinstance(value, str):
            self._suffix = value
        else:
            raise ValueError('Suffix must be a string.')
        # end if
    # end def

    def exists(self, filename: Path | str) -> bool:
        filename = self.get_filename(filename)
        return filename.exists()
    # end def

    def get_filename(self, filename: Path | str) -> Path:
        if isinstance(filename, str):
            filename = Path(filename)
        # end if
        if filename.is_file():
            filename = filename
        else:
            filename = filename / self.suffix
        # end if
        return filename
    # end def

    def load_result(
        self,
        filename: Path | str,
        default: list | ndarray | None | FileNotFoundError = FileNotFoundError(),
        rescale=True,
        **kwargs
    ) -> ndarray:
        filename = self.get_filename(filename)
        if filename.exists():
            data = loadtxt(filename, **kwargs)
            if data.ndim == 0:
                data = float(data)  # convert 0D array to float
            # end if
            # Normally rescale the data by the scale factor, but if rescale is False, return the raw data
            if rescale:
                data /= self.scale
            # end if
        else:
            if isinstance(default, FileNotFoundError):
                # If no default is given, raise an error if the file is missing
                raise FileNotFoundError(f"Could not find {filename}")
            else:
                data = default
            # end if
        # end if
        return data
    # end def

    def save_result(
        self,
        filename: Path | str,
        data: ndarray | float,
        overwrite=False,
        **kwargs
    ) -> None:
        filename = self.get_filename(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)
        if not filename.exists() or overwrite:
            if data is None:
                data = nan
            # end if
            if isscalar(data):
                data = array([data])
            # end if
            savetxt(filename, data, **kwargs)
        # end if
    # end def

# end class
