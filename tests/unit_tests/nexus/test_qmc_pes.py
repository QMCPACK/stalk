#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises

from stalk.nexus.qmc_pes import QmcPes
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_function import NotEvaluatedException
from stalk.util.util import match_to_tol


def test_QmcPes(tmp_path):

    # Test with empty args / defaults
    pes = QmcPes()
    assert isinstance(pes, QmcPes)
    assert pes.suffix == 'dmc/dmc.in.xml'
    assert pes.scale == 1.0

    # default suffix: dmc/dmc.in.xml
    # See tests/unit_tests/assets/qmc_pes/dmc/dmc.s001.scalar.dat
    E_ref = -37.630278
    Err_ref = 0.004255
    path = 'tests/unit_tests/assets/qmc_pes'
    s = ParameterSet()
    s.path = path
    pes.args = {}  # ensure that other tests are not interfering with the args
    res = pes.load(path)
    assert match_to_tol(res.value, E_ref, 1e-5)
    assert match_to_tol(res.error, Err_ref, 1e-5)

    # Test equilibration and qmc_idx=0
    # See tests/unit_tests/assets/qmc_pes/dmc/dmc.s000.scalar.dat
    # equilibration=25, ElecElec
    E1_ref = 116.856751
    Err1_ref = 0.018098
    pes.args = {}  # ensure that other tests are not interfering with the args
    pes.args['equilibration'] = 25
    pes.args['term'] = 'ElecElec'
    pes.args['qmc_idx'] = 0
    pes.scale = 3.0
    res1 = pes.load(path)
    assert match_to_tol(res1.value, E1_ref / 3, 1e-5)
    assert match_to_tol(res1.error, Err1_ref / 3, 1e-5)

    # Test failing analysis (e.g. missing energy term)
    pes.suffix = 'dmc_fail/dmc.in.xml'
    with raises(NotEvaluatedException):
        pes.load(path)
    # end with

    # Test skipping of missing file
    with raises(NotEvaluatedException):
        pes.load('missing')
    # end with

# end def
