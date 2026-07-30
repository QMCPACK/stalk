#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.nexus.pwscf_geometry import PwscfGeometry
from stalk.io.xyz_geometry import XyzGeometry
from stalk.params import GeometryResult
from stalk.params.pes_function import NotEvaluatedException
from stalk.util.util import match_to_tol


def test_PwscfGeometry():

    # Test with empty args
    pes = PwscfGeometry()
    assert pes.suffix == 'relax.in'
    assert pes.scale == 1.0

    # Use XyzLoader for reference
    pos_ref = XyzGeometry(suffix='relax_bohr.xyz').load(
        'tests/unit_tests/assets/pwscf_relax',
    )

    # default suffix: relax.in; only path is needed
    res0 = pes.load('tests/unit_tests/assets/pwscf_relax')
    assert isinstance(res0, GeometryResult)
    assert res0.get_axes() is None
    assert res0.get_elem() is None
    assert match_to_tol(res0.get_pos(), pos_ref.get_pos(), 1e-6)

    # Test by providing args (c_pos multiplication by 2)
    loader1 = PwscfGeometry(suffix='pwscf_relax/relax.in', c_pos=2.0)
    res1 = loader1.load('tests/unit_tests/assets')
    assert isinstance(res1, GeometryResult)
    assert match_to_tol(res1.get_pos(), 2.0 * pos_ref.get_pos(), 1e-6)

    # Test reading axes
    loader1.suffix = 'relax_axes.in'
    res2 = loader1.load('tests/unit_tests/assets/pwscf_relax')
    assert isinstance(res2, GeometryResult)
    # Only test superficially for now
    assert res2.get_axes() is not None

    # Test structures being missing (loading PES job)
    with raises(NotEvaluatedException):
        loader1.load('tests/unit_tests/assets/pwscf_pes/scf.in')
    # end with
# end def
