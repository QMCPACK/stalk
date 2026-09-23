#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from pathlib import Path
from numpy import ndarray

from stalk.io.txt_data import TxtData


class Cacheable():
    """Class for caching numerical attributes to disk using TxtData writers and readers."""
    _cache: dict[str, TxtData] = None

    def __init__(
        self,
        **data: dict[str, TxtData]
    ):
        self._cache = {}
        for key, value in data.items():
            self.add_cache(key, value)
        # end for
    # end def

    def add_cache(self, key: str, value: TxtData) -> None:
        if isinstance(value, TxtData):
            self._cache[key] = value
        else:
            raise TypeError(f"Value for key '{key}' must be a TxtData instance!")
        # end if
    # end def

    def save(self, path: str | Path | None = None, overwrite=False, **data) -> None:
        if path is None:
            return
        # end if
        for key, value in data.items():
            # Find the writer in the Cacheable configuration
            if key in self._cache.keys():
                writer = self._cache[key]
                if value is not None:
                    # The value must be float, list[float] or ndarray
                    writer.save_result(path, value, overwrite=overwrite)
                # end if
            else:
                # This should not happen if the inheriting classes are implemented right.
                raise KeyError(f"Field '{key}' is not cacheable!")
            # end if
        # end for
    # end def

    def load(
        self,
        path: str | Path | None = None,
        key: str = None,
        default=None
    ) -> dict[str, float | ndarray] | ndarray:
        """Load all or one specific cacheable data from the given path."""
        if path is None:
            return
        # end if
        if key is None:
            result = {}
            # Load all cacheable data
            for key in self._cache.keys():
                result[key] = self._cache[key].load_result(path, default=default)
            # end for
        else:
            if key in self._cache.keys():
                result = self._cache[key].load_result(path, default=default)
            else:
                raise KeyError(f"Field '{key}' is not in the cache!")
            # end if
        # end if
        return result
    # end def

# end class
