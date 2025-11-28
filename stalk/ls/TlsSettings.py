#!/usr/bin/env python3
'''Target line-search settings'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import isscalar, ndarray, random, where
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from stalk.ls.FittingResult import FittingResult
from stalk.ls.LsSettings import LsSettings


class TlsSettings(LsSettings):
    _Gs: ndarray = None
    _bias_mix: float
    _bias_order: int
    _target: FittingResult
    _interp = CubicHermiteSpline  # Interpolant

    def __init__(
        self,
        M=None,
        N=200,
        Gs=None,
        bias_order=1,
        bias_mix=0.0,
        target_x0=0.0,
        target_y0=0.0,
        **ls_args
        # fraction=0.025, sgn=1, fit_kind=None, fit_func=None, fit_args={},
    ):
        super().__init__(self, **ls_args)
        self.regenerate_Gs(M, N, Gs)
        self.bias_order = bias_order
        self.bias_mix = bias_mix
        self._target = FittingResult(target_x0, target_y0)
        self._interp = None
    # end def

    @property
    def bias_mix(self):
        return self._bias_mix
    # end def

    @bias_mix.setter
    def bias_mix(self, bias_mix):
        if isscalar(bias_mix) and bias_mix >= 0.0:
            self._bias_mix = bias_mix
        else:
            raise ValueError("Bias mix must be > 0, " + str(bias_mix) + " provided")
        # end if
    # end def

    @property
    def bias_order(self):
        return self._bias_order
    # end def

    @bias_order.setter
    def bias_order(self, bias_order):
        if isinstance(bias_order, int) and bias_order > 0:
            self._bias_order = bias_order
        else:
            raise ValueError("Bias order must be >0 integer")
        # end if
    # end def

    @property
    def M(self):
        if self.Gs is None:
            return 0
        else:
            return self.Gs.shape[1]
        # end if
    # end def

    @property
    def N(self):
        if self.Gs is None:
            return 0
        else:
            return self.Gs.shape[0]
        # end if
    # end def

    # Dummy setter is provided to override LsSettings with no effect
    @N.setter
    def N(self, N):
        pass
    # end def

    @property
    def Gs(self):
        return self._Gs
    # end def

    @property
    def interp(self):
        return self._interp
    # end def

    @interp.setter
    def interp(self, interp):
        self._interp = interp
    # end def

    @property
    def target(self):
        return self._target
    # end def

    @Gs.setter
    def Gs(self, Gs):
        if isinstance(Gs, ndarray) and Gs.ndim == 2:
            if Gs.shape[0] < 2:
                raise AssertionError('Gs must have at least two rows')
            elif Gs.shape[1] < 3:
                raise AssertionError('Gs must have at least three columns')
            else:
                self._Gs = Gs
            # end if
        elif Gs is None:
            self._Gs = None
        else:
            raise ValueError('Gs must be a 2D array')
        # end if
    # end def

    @property
    def interp_kind(self):
        if isinstance(self.interp, PchipInterpolator):
            return 'pchip'
        elif isinstance(self.interp, CubicHermiteSpline):
            return 'cubic'
        else:
            return None
        # end if
    # end def

    # Reset or regenerate Gs if the user changes M or N or provides new Gs
    def regenerate_Gs(self, M=None, N=None, Gs=None):
        if Gs is None and M is None and N is None:
            # Reset to None
            self.Gs = None
        elif Gs is not None:
            # If Gs provided, try to insert them and ignore other arguments
            self.Gs = Gs
        else:
            # Check if M or N change; if so, regenerate Gs accordingly or raise error
            M = M if M is not None else self.M
            N = N if N is not None else self.N
            if M != self.M and N != self.N:
                # The values of M and N are checked in Gs.setter
                self.Gs = random.randn(N, M)
            # end if
        # end if
    # end def

    def copy(
        self,
        Gs=None,
        M=None,
        N=None,
        **ls_overrides
    ):
        Gs = Gs if Gs is not None else self.Gs
        ls_args = {
            'fraction': self.fraction,
            'sgn': self.sgn,
            'fit_func': self.fit_func,
            'bias_mix': self.bias_mix,
            'bias_order': self.bias_order,
            'Gs': Gs,
        }
        M = M if M is not None else self.M
        N = N if N is not None else self.N
        # If no Gs/M/N are supplied, copy Gs from previous
        if M != self.M or N != self.N:
            ls_args['Gs'] = None
            ls_args['M'] = M
            if N > 0:
                ls_args['N'] = N
            # end if
        # end if
        ls_args.update(**ls_overrides)
        settings = TlsSettings(**ls_args)
        settings._interp = self.interp
        settings._target = self.target
        return settings
    # end def

    def get_safe_offsets(self, offsets):
        if self.interp is None or len(offsets) == 0:
            warnings.warn('Requested offsets but no interpolation is present!')
            return offsets
        # end if
        x_min = self.interp.x.min()
        x_max = self.interp.x.max()
        for i_offset in where(offsets < x_min)[0]:
            offset = offsets[i_offset]
            if self.interp.extrapolate:
                warnings.warn(f'Extrapolating for offset={offset} < R_min={x_min}')
            else:
                warnings.warn(f'Reset offset={offset} to R_min={x_min}')
                offsets[i_offset] = x_min
            # end if
        # end for
        for i_offset in where(offsets > x_max)[0]:
            offset = offsets[i_offset]
            if self.interp.extrapolate:
                warnings.warn(f'Extrapolating for offset={offset} < R_max={x_max}')
            else:
                warnings.warn(f'Reset offset={offset} to R_max={x_max}')
                offsets[i_offset] = x_max
            # end if
        # end for
        return offsets
    # end def

    # Return true if the settings of self and other are consistent
    def __eq__(self, other):
        if not isinstance(other, TlsSettings):
            return False
        # end if
        result = LsSettings.__eq__(self, other)
        result &= self.Gs is other.Gs
        result &= self.interp is other.interp
        result &= self.target is other.target
        result &= self.bias_order == other.bias_order
        result &= self.bias_mix == other.bias_mix
        return result
    # end def

# end class
