#!/usr/bin/env python3
'''PesResult represents a PES evaluation result as value+error pair (float/nan)'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import isscalar, nan, isnan

from stalk.util.noise import AbsNoise, Noise, WhiteNoise


class PesResult:
    _value = nan
    _error = 0.0

    @property
    def value(self):
        return self._value
    # end def

    @value.setter
    def value(self, value):
        if isscalar(value):
            self._value = value
        else:
            raise ValueError("Value must be scalar")
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
            raise ValueError("Error must be scalar and >= 0")
        # end if
    # end def

    def __init__(self, value, error=0.0):
        self.value = value
        self.error = error
    # end def

    def rescale(self, scale):
        self.value /= scale
        self.error /= scale
    # end def

    def add_sigma(self, sigma, kind: bool | Noise = True):
        '''Generate random noise to the result for error resampling purposes.'''
        if not kind:
            # False -> do not add noise
            return
        elif isinstance(kind, bool) and kind:
            # True -> use standard white noise
            noise = WhiteNoise()
        elif isinstance(kind, str):
            if kind.lower() == 'std':
                noise = WhiteNoise()
            elif kind.lower() == 'abs':
                noise = AbsNoise()
            # end if
        else:
            noise = kind
        # end if
        if not isinstance(noise, Noise):
            raise ValueError('Noise kind must be a Noise object or True/False')
        # end if
        if isinstance(sigma, float) and sigma >= 0.0:
            self.error = (self.error**2 + sigma**2)**0.5
            if not isnan(self.value):
                self.value += sigma * noise.generate_one()
            # end if
        else:
            raise ValueError('Tried to add poor sigma: ' + str(sigma))
        # end if
    # end def

# end class
