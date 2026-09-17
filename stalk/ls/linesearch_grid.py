#!/usr/bin/env python3
'''Class for containing a 1D grid of points, values and errorbars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from typing import TypeVar
import warnings
from matplotlib import pyplot as plt
from numpy import array, all, ndarray, searchsorted

from stalk.params.linesearch_point import LineSearchPoint
from stalk.pes.structure_collection import StructureCollection
from stalk.util.util import FFS

# Passing on the generic type
T = TypeVar('T', bound=LineSearchPoint)


class LineSearchGrid[T](StructureCollection[T]):

    def __init__(
        self,
        offsets=None,
        values=None,
        errors=None
    ):
        StructureCollection.__init__(self)
        if offsets is not None:
            if values is None:
                values = len(offsets) * [None]
                errors = len(offsets) * [0.0]
            elif errors is None:
                errors = len(offsets) * [0.0]
            # end if
            for offset, value, error in zip(offsets, values, errors):
                point = LineSearchPoint(offset, value, error)
                self.add_point(point)
            # end for
        # end for
    # end def

    @property
    def shifted(self) -> bool:
        '''True if more than two enabled points have been shifted'''
        return len([point for point in self._grid if point.enabled]) > 2
    # end def

    @property
    def evaluated(self) -> bool:
        '''True if all enabled points are evaluated'''
        return len(self) > 0 and all([point.valid for point in self._grid if point.enabled])
    # end def

    @property
    def valid_grid(self) -> ndarray:
        '''Return offset array of valid points'''
        return array([point for point in self._grid if point.valid])
    # end def

    @property
    def valid_args(self) -> ndarray:
        '''Return valid (offsets, values, errors) as a tuple of arrays'''
        return self.valid_offsets, self.valid_values, self.valid_errors
    # end def

    @property
    def offsets(self) -> ndarray:
        '''Return offset array of points'''
        return array([point.offset for point in self._grid])
    # end def

    @property
    def valid_offsets(self) -> ndarray:
        '''Return offset array of points'''
        return array([point.offset for point in self._grid if point.valid])
    # end def

    @property
    def valid_values(self) -> ndarray:
        '''Return values array of valid points'''
        return array([point.value for point in self._grid if point.valid])
    # end def

    @property
    def values(self) -> ndarray:
        '''Return values array of points'''
        return array([point.value for point in self._grid])
    # end def

    @values.setter
    def values(self, values) -> None:
        if len(values) == len(self):
            for value, point in zip(values, self._grid):
                point.value = value
            # end for
        else:
            raise ValueError("Values must be of same length as grid")
        # end if
    # end def

    @property
    def valid_errors(self) -> ndarray:
        '''Return errors array of valid points'''
        return array([point.error for point in self._grid if point.valid])
    # end def

    @property
    def errors(self) -> ndarray:
        '''Return errors array of valid points'''
        return array([point.error for point in self._grid])
    # end def

    @errors.setter
    def errors(self, errors) -> None:
        if len(errors) == len(self):
            for error, point in zip(errors, self._grid):
                point.error = error
            # end for
        else:
            raise ValueError("Errors must be of same length as grid")
        # end if
    # end def

    @property
    def R_max(self) -> float:
        if len(self) > 0:
            return min([-self.offsets.min(), self.offsets.max()])
        else:
            return 0.0
        # end if
    # end def

    @property
    def valid_R_max(self) -> float:
        if len(self.valid_grid) > 1:
            return min([-self.valid_offsets.min(), self.valid_offsets.max()])
        else:
            return 0.0
        # end if
    # end def

    @property
    def noisy(self) -> bool:
        '''True if any valid point has non-zero error'''
        return not all(self.valid_errors == 0.0)
    # end def

    @property
    def valid(self) -> bool:
        '''True if more than two valid points are present'''
        return len(self.valid_grid) > 2
    # end def

    def get(self, point: float | int | LineSearchPoint, default=None) -> LineSearchPoint | None:
        '''Get a requested point by the offset. Returns default if not found.'''
        # Find the first occurrence of the point in the grid, or return default if not found
        if isinstance(point, int) and abs(point) < len(self._grid):
            index = point
            return self._grid[index]
        # end if
        offset = point if isinstance(point, float) else point.offset
        index = searchsorted(self.offsets, offset)
        if index < len(self._grid) and self._grid[index].offset == offset:
            return self._grid[index]
        else:
            return default
        # end if
    # end def

    # Enable a point by the offset, if present
    def enable_value(self, offset) -> None:
        point = self.get(offset)
        if point is not None:
            point.enabled = True
        # end if
    # end def

    # Disable a point by offset, if present
    def disable_value(self, offset) -> None:
        point = self.get(offset)
        if point is not None:
            point.enabled = False
        # end if
    # end def

    def plot(
        self,
        ax=None,
        f=None,
        color='tab:blue',
        **kwargs
    ):
        if not self.valid:
            warnings.warn("Cannot plot without valid data.")
            return
        # end if
        if ax is None:
            f, ax = self._create_plot(**kwargs)
        # end if
        self.grid[0].plot(ax, color=color, label='Data', **kwargs)
        for point in self.grid[1:]:
            point.plot(ax, color=color, **kwargs)
        # end for
        plt.tight_layout()
    # end def

    def _create_plot(
        self,
        xlabel='Offset',
        ylabel='Energy',
        **kwargs
    ):
        f, ax = plt.subplots()
        ax.set_title(repr(self))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return f, ax
    # end def

    def __str__(self):
        string = repr(self)
        if len(self) == 0:
            string += '\nGrid: not set.'
        else:
            string += '\n  ' + (FFS + FFS + FFS).format('offset', 'value', 'error')
            for point in self.grid:
                string += '\n  ' + LineSearchPoint.__str__(point)
            # end for
        # end if
        return string
    # end def

    def __repr__(self):
        return self.__class__.__name__
    # end def

# end class
