#!/usr/bin/env python

from pytest import raises
from stalk.util import match_to_tol

from ..assets.h2o import pes_H2O, get_structure_H2O, get_hessian_H2O

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# Test ParallelLineSearch class
def test_ParallelLineSearch(tmp_path):
    from stalk import ParallelLineSearch

    with raises(TypeError):
        # Cannot init without PES
        ParallelLineSearch()
    # end with

    # Test empty init
    pls = ParallelLineSearch(pes_func=pes_H2O)
    assert not pls.setup
    assert not pls.shifted
    assert not pls.evaluated
    assert len(pls) == 0
    assert pls.D == 0
    assert pls.structure is None
    assert pls.hessian is None
    assert len(pls.ls_list) == 0
    assert len(pls.Lambdas) == 0
    assert len(pls.noises) == 0
    assert len(pls.windows) == 0
    assert not pls.noisy
    assert pls.params is None
    assert pls.params_err is None

    # Add structure
    s = get_structure_H2O()
    pls.structure = s
    assert not pls.setup
    # Add Hessian
    h = get_hessian_H2O()
    pls.hessian = h
    assert pls.setup
    assert not pls.shifted
    with raises(AssertionError):
        pls.evaluate()
    # end with
    # To shift structures and generate ls_list, provide windows, noises etc
    M = 5
    fit_kind = 'pf2'
    N = 100
    pls.initialize(
        window_frac=0.05,
        M=M,
        fit_kind=fit_kind,
        N=N
    )
    assert pls.shifted
    assert not pls.evaluated
    assert pls.structure_next is None

    pls.evaluate()
    assert pls.evaluated
    # Checking against hard-coded references
    assert match_to_tol(pls.structure_next.params, [1.02550264, 104.12792928])
    assert match_to_tol(pls.structure_next.params_err, [0.0, 0.0])

    # Test propagate and write
    pls.path = str(tmp_path)
    fname = 'test_fname.p'
    pls_next = pls.propagate(write=True, fname=fname)
    assert pls_next.path == str(tmp_path) + '_next/'

    # Test loading of pickle
    pls_load = ParallelLineSearch(load=str(tmp_path) + '/' + fname)
    assert pls_load.evaluated

# end def
