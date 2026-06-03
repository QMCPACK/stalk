#!/usr/bin/env python3

from pytest import raises

from stalk.nexus.pwscf_pes import PwscfPes
from stalk.params.pes_function import NotEvaluatedException
from stalk.util.util import match_to_tol


def test_PwscfPes(tmp_path):

    # Test with empty args / defaults
    pes = PwscfPes()
    assert pes.suffix == 'scf.in'
    assert pes.scale == 1.0

    # default suffix: scf.in
    E_ref = -22.74988263  # See tests/unit_tests/assets/pwscf_pes/scf.out
    pes.args = {}  # ensure that other tests are not interfering with the args
    res = pes.load('tests/unit_tests/assets/pwscf_pes')
    assert match_to_tol(res.value, E_ref)
    assert match_to_tol(res.error, 0.0)

    # failing output file
    with raises(NotEvaluatedException):
        pes.load('tests/unit_tests/assets/pwscf_pes/scf_failed.out')
    # end with

    # Test skipping of missing test
    with raises(NotEvaluatedException):
        pes.load('missing')
    # end with

# end def
