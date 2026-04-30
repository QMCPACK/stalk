#!/usr/bin/env python3
'''Class for containing grid, value and error of a single point of line-search'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from matplotlib.axes import Axes
from numpy import isscalar, abs, isnan

from stalk.util.util import FF, FFS, FU


class LineSearchPoint():
    _offset = 0.0
    _value = None
    _error = 0.0
    _enabled = True
    tol = 1e-9

    def __init__(
        self,
        offset,
        value=None,
        error=0.0
    ):
        self.offset = offset
        self.value = value
        self.error = error
    # end def

    # Offset relative to external reference identifies the point along grid
    @property
    def offset(self):
        return self._offset
    # end def

    @offset.setter
    def offset(self, offset):
        if not isscalar(offset):
            raise TypeError('Grid offset must be scalar')
        # end if
        self._offset = offset
    # end def

    @property
    def value(self):
        return self._value
    # end def

    @property
    def is_eqm(self):
        return abs(self.offset) < self.tol
    # end def

    @value.setter
    def value(self, value):
        if isscalar(value) or value is None:
            self._value = value
        else:
            raise ValueError('Value must be scalar or None')
        # end if
    # end def

    @property
    def error(self):
        return self._error
    # end def

    @error.setter
    def error(self, error):
        if isscalar(error) and error >= 0.0:
            self._error = error
        else:
            raise ValueError('Error must be a non-negative scalar.')
        # end if
    # end def

    @property
    def enabled(self):
        return self._enabled
    # end def

    @enabled.setter
    def enabled(self, enabled):
        self._enabled = enabled
    # end def

    @property
    def valid(self):
        '''The value is valid, when it is enabled and has a value'''
        return self._enabled and isscalar(self._value) and not isnan(self._value)
    # end def

    def reset_value(self):
        self.value = None
        self.error = 0.0
    # end def

    def plot(
        self,
        ax: Axes,
        color='tab:blue',
        marker='o',
        linestyle='none',
        **kwargs
    ):
        if self.error > 0:
            ax.errorbar(
                self.offset,
                self.value,
                self.error,
                color=color,
                marker=marker,
                linestyle=linestyle,
                **kwargs
            )
        else:
            ax.plot(
                self.offset,
                self.value,
                color=color,
                marker=marker,
                linestyle=linestyle,
                **kwargs
            )
        # end if
    # end def

    # Compare equality of two grid points
    def __eq__(self, other):
        return isinstance(other, LineSearchPoint) and abs(self.offset - other.offset) < self.tol
    # end def

    # Compare ordering of two grid points
    def __lt__(self, other):
        return isinstance(other, LineSearchPoint) and self.offset < other.offset
    # end def

    def __str__(self):
        string = FF.format(self.offset)
        if self.value is None:
            string += FFS.format('-')
        else:
            string += FF.format(self.value)
        # end if
        if self.error > 0:
            string += FU.format(self.error)
        else:
            string += FFS.format('-')
        # end if
        return string
    # end def

# end class
