#!/usr/bin/env python

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.params.pes_function import PesFunction
from stalk.pls.surrogate import Surrogate
from stalk.util import match_to_tol

from stalk.lsi import LineSearchIteration

from ..assets.h2o import get_structure_H2O, get_hessian_H2O, pes_H2O


# test LineSearchIteration class
def test_linesearchiteration(tmp_path):

    lsi = LineSearchIteration()
    assert len(lsi) == 0
    assert lsi.path == ''

    # Test default init from Hessian and structure
    pes = PesFunction(pes_H2O)
    path0 = str(tmp_path) + '/lsi0'
    hessian = get_hessian_H2O()
    structure = get_structure_H2O()
    structure.shift_params([0.2, -0.2])
    lsi = LineSearchIteration(
        path=path0,
        hessian=hessian,
        structure=structure,
    )
    assert len(lsi) == 1
    assert lsi.path == path0 + '/'
    # We can call evaluate explicitly
    lsi.evaluate(pes)
    assert lsi.evaluated
    # And then propagate
    # (defaults: fname='pls.p', write=True, overrite=True, add_sigma=False)
    lsi.propagate(pes)
    # Length should now be 2
    assert len(lsi) == 2
    assert lsi[0].evaluated
    assert not lsi[1].evaluated
    # We can readily propagate (evaluate is called therein)
    lsi.propagate(pes)
    assert len(lsi) == 3
    lsi.propagate(pes, add_sigma=True, write=False)
    assert len(lsi) == 4
    # Now, let's start by loading

    lsi_load = LineSearchIteration(path=path0)
    # The last iteration was not written
    assert len(lsi_load) == 2
    for i in range(len(lsi_load)):
        assert match_to_tol(lsi_load[i].structure.params, lsi[i].structure.params)
    # end for

    # Test default init from surrogate
    srg = Surrogate(
        fit_kind='pf3',
        hessian=hessian,
        structure=structure,
    )
    srg.evaluate(pes)
    with raises(AssertionError):
        # Cannot copy before optimized
        lsi_srg = LineSearchIteration(surrogate=srg)
    # end with
    windows = [0.1, 0.2]
    noises = [0.03, 0.04]
    M = 9
    srg.optimize_windows_noises(
        fit_kind='pf4',
        M=M,
        windows=windows,
        noises=noises
    )
    lsi_srg = LineSearchIteration(surrogate=srg)
    # Not the same object but same values
    assert lsi_srg[-1].structure is not srg.structure
    assert match_to_tol(lsi_srg[-1].structure.params, srg.structure.params)
    assert match_to_tol(lsi_srg[-1].hessian.hessian, srg.hessian.hessian)
    assert match_to_tol(lsi_srg[-1].windows, windows)
    assert match_to_tol(lsi_srg[-1].noises, noises)
    assert len(lsi_srg[-1].ls(0)) == M
    assert len(lsi_srg[-1].ls(1)) == M

# end def
