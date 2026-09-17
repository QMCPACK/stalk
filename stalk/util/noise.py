#!/usr/bin/env python3
"""PES noise model classes."""

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from numpy import random, ndarray

from stalk.util.args_container import ArgsContainer


class Noise(ArgsContainer):
    """Abstract class for a PES noise model"""

    def generate(self, N: int, M: int) -> ndarray:
        """Generate N by M noise value array"""
        raise NotImplementedError("Noise generation not implemented")
    # end def

    def generate_one(self) -> ndarray:
        """Generate noise"""
        raise NotImplementedError("Noise generation not implemented")
    # end def

# end class


class WhiteNoise(Noise):
    """White noise model"""

    def generate(self, N: int, M: int) -> ndarray:
        """Generate N by M standard white noise value array"""
        return random.normal(loc=0.0, scale=1.0, size=(N, M))
    # end def

    def generate_one(self) -> float:
        """Generate standard white noise"""
        return random.normal(loc=0.0, scale=1.0)
    # end def

# end class


class AbsNoise(WhiteNoise):
    """White noise model"""

    def generate(self, N: int, M: int) -> ndarray:
        """Generate N by M standard absolute white noise value array"""
        Gs = abs(super().generate(N=N, M=M))
        return Gs
    # end def

    def generate_one(self) -> float:
        """Generate standard absolute white noise"""
        return abs(super().generate_one())
    # end def

# end class


class NoiseFactory:
    """Factory class for creating Noise objects"""

    @staticmethod
    def create(kind: str | Noise) -> Noise:
        """Create a Noise object based on the kind"""
        if isinstance(kind, Noise):
            return kind
        elif kind == 'std':
            return WhiteNoise()
        elif kind == 'abs':
            return AbsNoise()
        else:
            raise ValueError(f"Unknown noise kind: {kind}")
        # end if
    # end def
# end class
