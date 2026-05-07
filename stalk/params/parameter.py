#!/usr/bin/env python
"""Base class for representing an optimizable parameters."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import isscalar, inf, ndarray, pi

from stalk.params.util import mean_bond_angles, mean_distances
from stalk.util.util import FF, FLL, FU


class Parameter():
    _value: float
    _error: float = 0.0
    _limits: tuple[float, float]
    label: str = ''
    unit: str = ''

    def __init__(
        self,
        value,
        error=0.0,
        label='p',
        unit='',
        limits=(-inf, inf),
    ):
        self.limits = limits
        self.value = value
        self.error = error
        self.label = label
        self.unit = unit
    # end def

    @property
    def value(self):
        return self._value
    # end def

    @value.setter
    def value(self, value):
        if isscalar(value):
            if value > self.limits[1]:
                raise ParameterLimitException(f"Cannot raise {self.label} above {self.limits[1]}")
            # end if
            if value < self.limits[0]:
                raise ParameterLimitException(f"Cannot decrease {self.label} below {self.limits[0]}")
            # end if
            self._value = value
        else:
            raise TypeError("Value must be scalar!")
        # end if
    # end def

    @property
    def error(self):
        return self._error
    # end def

    @error.setter
    def error(self, error):
        if isscalar(error):
            self._error = error
        else:
            self._error = 0.0
        # end if
    # end def

    @property
    def limits(self):
        return self._limits
    # end def

    @limits.setter
    def limits(self, limits):
        if limits is None:
            self._limits = (-inf, inf)
        elif isinstance(limits, (tuple, list)) and len(limits) == 2:
            self._limits = tuple(limits)
        else:
            raise TypeError(f'Limits must be tuple[float, float], provided: {limits}')
        # end if
    # end def

    def shift(self, shift: float):
        self.value += shift
        self.error = 0.0
    # end def

    def __str__(self):
        string = (FLL + FF).format(self.label, self.value)
        if self.error > 0:
            string += FU.format(self.error)
        # end if
        if self.unit is not None:
            string += FLL.format(self.unit)
        # end if
        return string
    # end def
# end class


class BondLength(Parameter):

    def __init__(
        self,
        pos: tuple[ndarray, ndarray] | list[tuple[ndarray, ndarray]],
        label='d',
        unit='A',
        error=0.0,
        limits=(0.0, inf),
        tol=1e-6,
        axes=None,
    ):
        if isinstance(pos, tuple):
            # Provide just one pair
            d = mean_distances([pos], tol=tol, axes=axes)
        else:
            # Provide list of pairs
            d = mean_distances(pos, tol=tol, axes=axes)
        # end if
        Parameter.__init__(
            self,
            value=d,
            error=error,
            label=label,
            unit=unit,
            limits=limits,
        )
    # end def

# end def


class BondAngle(Parameter):

    def __init__(
        self,
        pos: tuple[ndarray, ndarray, ndarray] | list[tuple[ndarray, ndarray, ndarray]],
        label='a',
        unit='ang',
        error=0.0,
        limits=(0.0, 180.0),
        tol=1e-6,
        axes=None,
    ):
        if isinstance(pos, tuple):
            # Provide just one pair
            a = mean_bond_angles([pos], tol=tol, axes=axes, units=unit)
        else:
            # Provide list of pairs
            a = mean_bond_angles(pos, tol=tol, axes=axes, units=unit)
        # end if
        Parameter.__init__(
            self,
            value=a,
            error=error,
            label=label,
            unit=unit,
            limits=limits,
        )
    # end def

# end def


class PhaseAngle(Parameter):

    def __init__(
        self,
        value,
        label='t',
        unit='rad',
        error=0.0,
    ):
        if unit == 'rad':
            limits = (-pi, pi)
        else:
            limits = (-180.0, 180.0)
        # end if
        Parameter.__init__(
            self,
            value=value,
            error=error,
            label=label,
            unit=unit,
            limits=limits,
        )
    # end def

    @Parameter.value.setter
    def value(self, value):
        if isscalar(value):
            # Reset between limits (works in both units as long as limits are proper)
            dlim = self.limits[1] - self.limits[0]
            value = (value + self.limits[1]) % dlim + self.limits[0]
            self._value = value
        else:
            raise TypeError("Value must be scalar!")
        # end if
    # end def

# end def


class ParameterLimitException(Exception):

    def __init__(self, msg):
        super().__init__(self, msg)
    # end def

# end class
