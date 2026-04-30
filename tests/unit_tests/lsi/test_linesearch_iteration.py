#!/usr/bin/env python

from pytest import raises
from stalk.pls.surrogate import Surrogate
from stalk.util import match_to_tol

from stalk.lsi import LineSearchIteration

from ..assets.h2o import get_structure_H2O, get_hessian_H2O, pes_H2O

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# test LineSearchIteration class
def test_linesearchiteration(tmp_path):

    # Test empty init
    with raises(TypeError):
        lsi = LineSearchIteration()
    # end with
    lsi = LineSearchIteration(pes_func=pes_H2O)
    assert len(lsi) == 0
    assert lsi.path == ''

    # Test default init from Hessian and structure
    path0 = str(tmp_path) + '/lsi0'
    hessian = get_hessian_H2O()
    structure = get_structure_H2O()
    structure.shift_params([0.2, -0.2])
    lsi = LineSearchIteration(
        path=path0,
        hessian=hessian,
        structure=structure,
        pes_func=pes_H2O,
    )
    assert len(lsi) == 1
    assert lsi.path == path0 + '/'
    # We can call evaluate explicitly
    lsi.evaluate()
    assert lsi.pls().evaluated
    # And then propagate
    # (defaults: fname='pls.p', write=True, overrite=True, add_sigma=False)
    lsi.propagate()
    # Length should now be 2
    assert len(lsi) == 2
    assert lsi.pls(0).evaluated
    assert not lsi.pls(1).evaluated
    # We can readily propagate (evaluate is called therein)
    lsi.propagate()
    assert len(lsi) == 3
    lsi.propagate(add_sigma=True, write=False)
    assert len(lsi) == 4
    # Now, let's start by loading

    lsi_load = LineSearchIteration(
        path=path0
    )
    # The last iteration was not written
    assert len(lsi_load) == 2
    for i in range(len(lsi_load)):
        assert match_to_tol(lsi_load.pls(i).structure.params, lsi.pls(i).structure.params)
    # end for

    # Test default init from surrogate
    srg = Surrogate(
        fit_kind='pf3',
        hessian=hessian,
        structure=structure,
        pes_func=pes_H2O
    )
    with raises(AssertionError):
        # Cannot copy before optimized
        lsi_srg = LineSearchIteration(surrogate=srg, pes_func=pes_H2O)
    # end with
    windows = [0.1, 0.2]
    noises = [0.03, 0.04]
    srg.optimize_windows_noises(
        fit_kind='pf4',
        M=7,
        windows=windows,
        noises=noises
    )
    lsi_srg = LineSearchIteration(surrogate=srg, pes_func=pes_H2O)
    # Not the same object but same values
    assert lsi_srg.pls().structure is not srg.structure
    assert match_to_tol(lsi_srg.pls().structure.params, srg.structure.params)
    assert match_to_tol(lsi_srg.pls().hessian.hessian, srg.hessian.hessian)
    assert match_to_tol(lsi_srg.pls().windows, windows)
    assert match_to_tol(lsi_srg.pls().noises, noises)
    assert len(lsi_srg.pls().ls(0)) == 7
    assert len(lsi_srg.pls().ls(1)) == 7

# end def
