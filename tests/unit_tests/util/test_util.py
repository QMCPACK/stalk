#!/usr/bin/env python

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array, ones, random
from pytest import raises

from stalk.util import match_to_tol, get_fraction_error


def test_units():
    from stalk.util.util import Bohr, Ry, Hartree
    assert match_to_tol(Bohr, 0.5291772105638411, 1e-10)
    assert match_to_tol(Ry, 13.605693012183622, 1e-10)
    assert match_to_tol(Hartree, 2 * Ry, 1e-10)
# end def


def test_get_fraction_error():

    # Test data
    N = array(range(101)) + 50
    N_skew = array(range(101)) + 50
    # Making one side '99's does not affect the median nor the errorbar
    N_skew[0:50] = 99 * ones(50)
    random.shuffle(N)
    random.shuffle(N_skew)

    # Test degraded
    with raises(ValueError):
        err, ave_b = get_fraction_error(N, fraction=0.5)
    # end with
    with raises(ValueError):
        err, ave_b = get_fraction_error(N, fraction=-1e-10)
    # end with

    # Test nominal
    for frac, target in zip([0.0, 0.01, 0.1, 0.49], [50, 49, 40, 1]):
        ave_b, err_b = get_fraction_error(N_skew, fraction=frac, both=True)
        ave, err = get_fraction_error(N_skew, fraction=frac, both=False)
        assert match_to_tol(err, target)
        assert match_to_tol(ave, 100)
        assert match_to_tol(err_b[0], 1)
        assert match_to_tol(err_b[1], target)
        assert match_to_tol(ave_b, 100)
    # end for

# end def


def test_match_to_tol():
    assert match_to_tol([[0.0, -0.1]], [0.1, -0.2], 0.1 + 1e-8)
# end def
