#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.io.geometry_loader import GeometryLoader


def test_GeometryLoader():

    gl = GeometryLoader()
    assert gl.scale == 1.0
    assert gl.suffix == 'structure.dat'

    gl = GeometryLoader(suffix='geom.dat', scale=0.5)
    assert gl.scale == 0.5

    gl = GeometryLoader(suffix='geom.dat', scale=0.5, c_pos=0.5)
    assert gl.scale == 2.0

    with raises(NotImplementedError):
        # File exists but its reading is not implemented
        gl.load('tests/unit_tests/assets/pwscf_relax/relax_bohr.xyz')
    # end with

# end def
