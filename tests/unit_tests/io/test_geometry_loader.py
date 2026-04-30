#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.io.geometry_loader import GeometryLoader


def test_GeometryLoader():

    gl = GeometryLoader()
    with raises(FileNotFoundError):
        gl.load('missing')
    # end with

    with raises(NotImplementedError):
        # File exists but its readin is not implemented
        gl.load('tests/unit_tests/assets/pwscf_relax/relax_bohr.xyz')
    # end with

# end def
