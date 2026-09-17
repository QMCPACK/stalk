#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.util.noise import Noise
from stalk.util.noise import WhiteNoise
from stalk.util.noise import AbsNoise
from stalk.util.noise import NoiseFactory


def test_Noise():

    noise = Noise()
    with raises(NotImplementedError):
        noise.generate(N=10, M=5)
    # end with

    with raises(NotImplementedError):
        noise.generate_one()
    # end with

# end def


def test_WhiteNoise():

    noise = WhiteNoise()
    Gs = noise.generate(N=10, M=5)
    assert Gs.shape == (10, 5)
    assert isinstance(noise.generate_one(), float)

# end def


def test_AbsNoise():

    noise = AbsNoise()
    Gs = noise.generate(N=10, M=5)
    assert Gs.shape == (10, 5)
    assert (Gs >= 0).all()
    assert isinstance(noise.generate_one(), float)

# end def


def test_NoiseFactory():

    noise = NoiseFactory.create('std')
    assert isinstance(noise, WhiteNoise)

    noise = NoiseFactory.create('abs')
    assert isinstance(noise, AbsNoise)

    with raises(ValueError):
        noise = NoiseFactory.create('invalid')
    # end with

    noise = Noise()
    assert noise is NoiseFactory.create(noise)

# end def
