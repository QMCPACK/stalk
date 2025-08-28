#!/usr/bin/env python3
"""Class for PES line-search in structure context"""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from matplotlib import pyplot as plt
from numpy import array, polyval, sign, isscalar

from stalk.ls.FittingResult import FittingResult
from stalk.params.ParameterHessian import ParameterHessian
from stalk.params.PesFunction import PesFunction
from stalk.params.ParameterSet import ParameterSet
from stalk.ls.LineSearchBase import LineSearchBase
from stalk.util.util import FF, SL


class LineSearch(LineSearchBase):
    _structure: ParameterSet = None  # The equilibrium structure
    _hessian: ParameterHessian = None  # The equilibrium full Hessian
    _sigma = 0.0  # Target errorbar
    _d: int = None  # direction count
    _enabled = True  # whether enabled or not

    def __init__(
        self,
        structure=None,
        hessian=None,
        d=None,
        sigma=0.0,
        offsets=None,
        M=7,
        W=None,
        R=None,
        pes=None,
        path='',
        dep_jobs=None,
        interactive=False,
        **ls_args
        # values=None, errors=None, fraction=0.025, sgn=1
        # fit_kind='pf3', fit_func=None, fit_args={}, N=200, Gs=None
    ):
        LineSearchBase.__init__(self, offsets=None, **ls_args)
        self.sigma = sigma
        if d is not None:
            self.d = d
        # end if
        if structure is not None:
            self.structure = structure
        # end if
        if hessian is not None:
            self.hessian = hessian
        # end if
        # Try to initialize grid based on available information
        try:
            self.set_grid(M=M, W=W, R=R, offsets=offsets)
            # Try to evaluate the pes and set the results
            if isinstance(pes, PesFunction):
                self.evaluate(
                    pes=pes,
                    interactive=interactive,
                    path=path,
                    dep_jobs=dep_jobs
                )
            # end if
        except (ValueError):
            # If the grid or pes input values are missing, the grid will be set later
            pass
        # end try
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
    def sigma(self):
        return self._sigma
    # end def

    @sigma.setter
    def sigma(self, sigma):
        if isscalar(sigma) and sigma >= 0.0:
            self._sigma = sigma
        else:
            raise ValueError("Sigma must be >= 0.0")
        # end if
    # end def

    @property
    def d(self):
        return self._d
    # end def

    @d.setter
    def d(self, d):
        if self.hessian is not None:
            max_d = len(self.hessian)
        elif self.structure is not None:
            max_d = len(self.structure)
        else:
            max_d = 1e10
        # end if
        if isinstance(d, int) and d < max_d:
            self._d = d
        else:
            raise ValueError('d must be integer smaller than the Hessian/structure dimension')
        # end if
    # end def

    @property
    def hessian(self):
        return self._hessian
    # end def

    @hessian.setter
    def hessian(self, hessian):
        if isinstance(hessian, ParameterHessian):
            self._hessian = hessian
            Lambda = self.hessian.lambdas[self.d]
            self.sgn = int(sign(Lambda))
            if self.structure is None:
                # Use Hessian structure if no other has been provided yet
                self.structure = hessian.structure
            # end if
        else:
            raise ValueError('Provided Hessian is not a ParameterHessian object')
        # end if
    # end def

    @property
    def Lambda(self):
        return None if self.hessian is None else abs(self.hessian.lambdas[self.d])
    # end def

    @property
    def direction(self):
        if self.d is None:
            return 0.0
        # end if
        if self.hessian is not None:
            return self.hessian.directions[self.d]
        elif self.structure is not None:
            # Get pure parameter direction
            direction = len(self.structure) * [0.0]
            direction[self.d] += 1.0
            return array(direction)
        else:
            return 0.0
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
    def W_max(self):
        return self._R_to_W(self.R_max)
    # end def

    @property
    def valid_W_max(self):
        return self._R_to_W(self.valid_R_max)
    # end def

    def figure_out_offsets(self, M=7, W=None, R=None, offsets=None):
        if offsets is not None:
            return offsets
        # end
        if M < 0:
            raise ValueError("Grid size M must be positive!")
        # end if
        if R is not None:
            return self._make_offsets_R(R, M=M)
        elif W is not None and self.hessian is not None:
            if self.hessian is None:
                raise ValueError('Must set Hessian before using W to set grid.')
            else:
                return self._make_offsets_W(W, M=M)
            # end if
        else:
            raise ValueError('Must provide grid, R or W to characterize grid.')
        # end if
    # end def

    def set_grid(
        self,
        **grid_kwargs  # M=7, W=None, R=None, offsets=None
    ):
        offsets = self.figure_out_offsets(**grid_kwargs)
        if self.structure is None:
            # Reverting to LineSearchPoints
            self.grid = offsets
        else:
            # Using shifted parametric structures
            self.grid = [self._shift_structure(offset) for offset in offsets]
        # end if
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
        if self.Lambda is None:
            return None
        else:
            return 0.5 * self.Lambda * R**2
        # end if
    # end def

    def add_shift(self, shift):
        if self.structure is None:
            # Reverting to LineSearchPoint
            self.add_point(shift)
        else:
            structure = self._shift_structure(shift)
            self.add_point(structure)
        # end if
    # end def

    def _shift_structure(self, shift):
        structure = self.structure.copy(offset=shift)
        if structure.is_eqm:
            # i.e. abs(offset) < threshold
            structure.label = 'eqm'
        else:
            label = SL.format(self.d, shift)
            structure.label = label
            structure.shift_params(shift * self.direction)
        # end if
        return structure
    # end def

    def evaluate(
        self,
        pes: PesFunction = None,
        add_sigma=False,
        **kwargs,  # path='', interactive=False, dep_jobs=None
    ):
        '''Evaluate the PES on the line-search grid using an evaluation function.'''
        pes.evaluate_all(
            self._grid,
            sigmas=len(self) * [self.sigma],
            add_sigma=add_sigma,
            **kwargs
        )
        self._search_and_store()
    # end def

    @property
    def shifted_params(self):
        if len(self) > 0:
            return array([structure.params for structure in self.grid if isinstance(structure, ParameterSet)])
        else:
            return None
        # end if
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
            a = 0.5 * self.sgn * self.Lambda
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
        return f'#{self.d} {super().__repr__()}'
    # end def

# end class
