#!/usr/bin/env python3

from numpy import isnan
from pytest import warns

from stalk.nexus.QmcPes import QmcPes
from stalk.params.ParameterSet import ParameterSet
from stalk.util.util import match_to_tol


def test_QmcPes(tmp_path):

    # Test with empty args / defaults
    pes = QmcPes()
    assert pes.func is None
    assert len(pes.args) == 0

    # default suffix: dmc/dmc.in.xml
    # See tests/unit_tests/assets/qmc_pes/dmc/dmc.s001.scalar.dat
    E_ref = -37.630278
    Err_ref = 0.004255
    s = ParameterSet()
    s.file_path = 'tests/unit_tests/assets/qmc_pes'
    res = pes.load(s)
    assert match_to_tol(res.value, E_ref, 1e-5)
    assert match_to_tol(res.error, Err_ref, 1e-5)

    # Test equilibration and qmc_idx=0
    # See tests/unit_tests/assets/qmc_pes/dmc/dmc.s000.scalar.dat
    # equilibration=25, ElecElec
    E1_ref = 116.856751
    Err1_ref = 0.018098
    res1 = pes.load(
        s,
        qmc_idx=0,
        equilibration=25,
        term='ElecElec'
    )
    assert match_to_tol(res1.value, E1_ref, 1e-5)
    assert match_to_tol(res1.error, Err1_ref, 1e-5)

    # failing output file
    with warns(UserWarning):
        res2 = pes.load(s, suffix='dmc_fail/dmc.in.xml')
        assert isnan(res2.value)
        assert match_to_tol(res2.error, 0.0)
    # end with

    # Test skipping of missing test
    with warns(UserWarning):
        s.file_path = 'missing'
        res_missing = pes.load(s)
        assert isnan(res_missing.value)
        assert match_to_tol(res_missing.error, 0.0)
    # end with

# end def
