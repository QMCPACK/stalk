#!/usr/bin/env python3
'''Class for containing a 1D grid of points, values and errorbars'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from matplotlib import pyplot as plt
from numpy import array, all

from stalk.params.linesearch_point import LineSearchPoint
from stalk.util.util import FFS


class LineSearchGrid():
    # List of LineSearchPoint instances
    _grid: list[LineSearchPoint] = []

    def __init__(
        self,
        offsets=None,
        values=None,
        errors=None
    ):
        self._grid = []

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
    def shifted(self):
        '''True if more than two enabled points have been shifted'''
        return len([point for point in self._grid if point.enabled]) > 2
    # end def

    @property
    def evaluated(self):
        '''True if all enabled points are evaluated'''
        return len(self) > 0 and all([point.valid for point in self._grid if point.enabled])
    # end def

    @property
    def valid_grid(self):
        '''Return offset array of valid points'''
        return array([point for point in self._grid if point.valid])
    # end def

    @property
    def grid(self):
        '''Return list of points'''
        return [point for point in self._grid]
    # end def

    @grid.setter
    def grid(self, grid):
        self._grid = []
        for point in grid:
            self.add_point(point)
        # end for
    # end def

    @property
    def offsets(self):
        '''Return offset array of points'''
        return array([point.offset for point in self._grid])
    # end def

    @property
    def valid_offsets(self):
        '''Return offset array of points'''
        return array([point.offset for point in self._grid if point.valid])
    # end def

    @property
    def valid_values(self):
        '''Return values array of valid points'''
        return array([point.value for point in self._grid if point.valid])
    # end def

    @property
    def values(self):
        '''Return values array of points'''
        return array([point.value for point in self._grid])
    # end def

    @values.setter
    def values(self, values):
        if len(values) == len(self):
            for value, point in zip(values, self._grid):
                point.value = value
            # end for
        else:
            raise ValueError("Values must be of same length as grid")
        # end if
    # end def

    @property
    def valid_errors(self):
        '''Return errors array of valid points'''
        return array([point.error for point in self._grid if point.valid])
    # end def

    @property
    def errors(self):
        '''Return errors array of valid points'''
        return array([point.error for point in self._grid])
    # end def

    @errors.setter
    def errors(self, errors):
        if len(errors) == len(self):
            for error, point in zip(errors, self._grid):
                point.error = error
            # end for
        else:
            raise ValueError("Errors must be of same length as grid")
        # end if
    # end def

    @property
    def R_max(self):
        if len(self) > 0:
            return min([-self.offsets.min(), self.offsets.max()])
        else:
            return 0.0
        # end if
    # end def

    @property
    def valid_R_max(self):
        if len(self) > 0 and len(self.valid_grid) > 1:
            return min([-self.valid_offsets.min(), self.valid_offsets.max()])
        else:
            return 0.0
        # end if
    # end def

    @property
    def noisy(self):
        return not all(self.valid_errors == 0.0)
    # end def

    @property
    def valid(self):
        return len(self.valid_grid) > 1
    # end def

    def add_point(self, point):
        if not isinstance(point, LineSearchPoint):
            point = LineSearchPoint(point)
        # end if
        if self.find_point(point) is None:
            self._grid.append(point)
            # Keep the grid sorted
            self._grid.sort()
        # end if
    # end def

    # Get all enabled arrays: (offsets, values, errors)
    def get_all(self):
        return self.offsets, self.values, self.errors
    # end def

    # Get full arrays including disabled and invalid: (offsets, values, errors)
    def get_valid(self):
        return self.valid_offsets, self.valid_values, self.valid_errors
    # end def

    # Finds and returns a requested point, if present; if not, returns None
    def find_point(self, point):
        # Try to find point by index
        if isinstance(point, int):
            if abs(point) < len(self):
                point = self._grid[point]
            else:
                return None
            # end if
        elif not isinstance(point, LineSearchPoint):
            # point is assumed to be a scalar offset
            point = LineSearchPoint(point)
        # end if
        # Else: point must be a LineSearchPoint
        for point_this in self._grid:
            if point_this == point:
                return point_this
            # end if
        # end for
    # end def

    # Sets the value and error for a given point, if present
    def set_value_error(self, offset, value, error=0.0):
        point = self.find_point(offset)
        if point is not None:
            point.value = value
            point.error = error
        # end if
    # end def

    # Enable a point by offset, if present
    def enable_value(self, offset):
        point = self.find_point(offset)
        if point is not None:
            point.enabled = True
        # end if
    # end def

    # Disable a point by offset, if present
    def disable_value(self, offset):
        point = self.find_point(offset)
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

    def __len__(self):
        return len(self.grid)
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
