#!/usr/bin/env python3
"""Class for PES line-search in structure context"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from matplotlib import pyplot as plt
from numpy import array, polyval, ndarray, zeros

from stalk.fit.fitting_result import FittingResult
from stalk.params.parameter_hessian import ParameterHessian
from stalk.params.parameter_set import ParameterSet
from stalk.ls.linesearch_base import LineSearchBase
from stalk.util.util import FF


class LineSearch(LineSearchBase):
    _structure: ParameterSet = None  # The equilibrium structure
    _direction: ndarray = None  # The search direction
    _Lambda: float = None  # The search direction stiffness
    _d: int = None  # Direction index

    def __init__(
        self,
        structure: ParameterSet = None,
        direction: ndarray = None,
        hessian: ParameterHessian = None,
        Lambda: float = None,
        d=None,
        offsets=None,
        values=None,
        errors=None,
        M=7,
        R=None,
        W=None,
        **ls_args
        # fraction=0.025, sgn=1, sigma=0
        # fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None
    ):
        # Hessian can be used to override structure, direction and Lambda
        if isinstance(hessian, ParameterHessian) and d is not None:
            structure = hessian.structure
            Lambda = hessian.lambdas[d]
            direction = hessian.directions[d]
            d = d
        # end if
        if structure is not None:
            self.structure = structure
            # Only set direction if the structure is set
            if direction is not None:
                self.direction = array(direction)
            elif d is not None:
                # Allow the direction to be a parameter direction
                direction = zeros(len(structure))
                direction[d] = 1.0
                self.direction = direction
            # end if
        # end if
        self.Lambda = Lambda
        self.d = d
        if offsets is None:
            if values is not None:
                warnings.warn("Grid offsets are automatically generated, ignoring provided values.")
                values = None
            # end if
            if errors is not None:
                warnings.warn("Grid offsets are automatically generated, ignoring provided errors.")
                errors = None
            # end if
            # Parent init ignoring offsets/values/errors
            LineSearchBase.__init__(self, **ls_args)
            # Finally, reset offsets
            self.grid = self.figure_out_offsets(M=M, W=W, R=R)
        else:
            # Default parent init
            LineSearchBase.__init__(
                self,
                offsets=offsets,
                values=values,
                errors=errors,
                **ls_args
            )
        # end if
    # end def

    @property
    def structure(self):
        return self._structure
    # end def

    @structure.setter
    def structure(self, structure):
        if isinstance(structure, ParameterSet):
            if structure.check_consistency():
                self._structure = structure
                # Empty grid when updating structure
                self._grid = []
            else:
                raise ValueError('Provided structure is not a consistent mapping')
            # end if
        else:
            raise ValueError('Provided structure is not a ParameterSet object')
        # end if
    # end def

    @property
    def Lambda(self):
        return self._Lambda
    # end def

    @Lambda.setter
    def Lambda(self, Lambda: float | None) -> None:
        if Lambda is not None:
            if Lambda <= 0.0:
                raise ValueError('Lambda must be positive')
            # end if
            self._Lambda = Lambda
        else:
            self._Lambda = None
        # end if
    # end def

    @property
    def direction(self):
        return self._direction
    # end def

    @direction.setter
    def direction(self, direction: ndarray) -> None:
        if self.structure is None:
            raise ValueError('Cannot set direction without a structure')
        # end if
        if isinstance(direction, ndarray):
            if len(direction) == len(self.structure):
                # Assume the direction is already normalized
                self._direction = direction
            else:
                raise ValueError('Direction must be of same length as structure')
            # end if
        else:
            raise TypeError('Direction must be a numpy array')
        # end if
    # end def

    @property
    def d(self) -> int | None:
        return self._d
    # end def

    @d.setter
    def d(self, d: int | None) -> None:
        if d is None or isinstance(d, int):
            self._d = d
        else:
            raise ValueError('d must be an integer or None')
        # end if
    # end def

    # Override to handle directional shifts
    @LineSearchBase.grid.setter
    def grid(self, grid: list[float]) -> None:
        self._grid = []
        for point in grid:
            if isinstance(point, float):
                self.add_shift(point)
            else:
                raise ValueError("Grid must be a list of floats")
            # end if
        # end for
    # end def

    @property
    def W_max(self):
        return self._R_to_W(self.R_max)
    # end def

    @property
    def valid_W_max(self):
        return self._R_to_W(self.valid_R_max)
    # end def

    @property
    def shifted_params(self):
        if len(self) > 0:
            return array([structure.params for structure in self.grid if isinstance(structure, ParameterSet)])
        else:
            return None
        # end if
    # end def

    def reset_offsets(self, M=7, W=None, R=None) -> None:
        self.grid = self.figure_out_offsets(M=M, W=W, R=R)
    # end def

    def figure_out_offsets(self, M=7, W=None, R=None) -> ndarray:
        if M < 0:
            raise ValueError("Grid size M must be positive!")
        # end if
        if R is not None:
            offsets = self._make_offsets_R(R, M=M)
        elif self.Lambda is None or W is None:
            offsets = []
        else:
            offsets = self._make_offsets_W(W, M=M)
        # end if
        if len(offsets) > 0 and self.structure is None:
            print('Could not reset offsets because structure and direction are not set.')
            offsets = []
        # end if
        return offsets
    # end def

    def add_shift(self, shift: float) -> None:
        if self.structure is None or self.direction is None:
            raise AssertionError("Cannot add shift without structure and direction set.")
        # end if
        structure = self._shift_structure(shift)
        self.add_point(structure)
    # end def

    def _make_offsets_W(self, W, M):
        if W < 0:
            raise ValueError("W must be positive!")
        # end if
        R = self._W_to_R(max(W, 1e-8))
        return self._make_offsets_R(R, M=M)
    # end def

    def _W_to_R(self, W):
        """Map W to R"""
        if self.Lambda is None:
            return None
        else:
            return (2 * W / self.Lambda)**0.5
        # end if
    # end def

    def _R_to_W(self, R):
        """Map R to W"""
        if self.Lambda is None or R is None:
            return None
        else:
            return 0.5 * self.Lambda * R**2
        # end if
    # end def

    def _shift_structure(self, shift):
        structure = self.structure.copy(offset=shift)
        if structure.is_eqm:
            # i.e. abs(offset) < threshold
            structure.label = 'eqm'
        else:
            if self.d is not None:
                structure.label = f'd{self.d}_{shift:+5.4f}'
            # end if
            structure.shift_params(shift * self.direction)
        # end if
        return structure
    # end def

    def plot(
        self,
        ax=None,
        color='tab:blue',
        target=None,
        **kwargs
    ):
        if not self.valid:
            warnings.warn("Cannot plot without valid data.")
            return
        # end if
        if ax is None:
            f, ax = self._create_plot(**kwargs)
        # end if
        if target is None:
            if self.fit_res is None:
                target = FittingResult(0.0, 0.0)
            else:
                target = self.fit_res
            # end if
        # end if
        LineSearchBase.plot(self, ax=ax, target=target, color=color, **kwargs)
        if self.Lambda is not None:
            a = 0.5 * self.settings.sgn * self.Lambda
            x0 = target.x0
            y0 = target.y0
            pfl = [a, -2 * a * x0, y0 + a * x0**2]
            xgrid = self._get_plot_grid(0.0)
            ygrid = polyval(pfl, xgrid)
            ax.plot(
                xgrid,
                ygrid,
                color=color,
                linestyle=':',
                label='Hessian'
            )
        # end if
        plt.tight_layout()
    # end def

    def __str__(self):
        string = LineSearchBase.__str__(self)
        if self.Lambda is not None:
            string += ('\n  Lambda: ' + FF).format(self.Lambda)
        # end if
        return string
    # end def

    def __repr__(self):
        if self.d is None:
            return super().__repr__()
        else:
            return f'#{self.d} {super().__repr__()}'
        # end if
    # end def

# end class
