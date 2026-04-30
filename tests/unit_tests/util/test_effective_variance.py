#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from numpy import isnan, array
from stalk.params.effective_variance import EffectiveVariance
from stalk.util.util import match_to_tol


def test_EffectiveVariance():

    ev0 = EffectiveVariance()
    with raises(ValueError):
        # Cannot add zero samples
        ev0.add_errorbar_data(0, 0.1)
    # end with
    with raises(ValueError):
        # Cannot add zero errorbar
        ev0.add_errorbar_data(0.1, 0)
    # end with
    assert isnan(ev0.var_eff)
    assert ev0.get_samples(2.0) == -1
    assert isnan(ev0.get_errorbar(2))

    # Test basic properties using scalar init
    samples = 5
    errorbar = 6.0
    ev1 = EffectiveVariance(samples, errorbar)
    # Test that we recover original samples
    assert ev1.get_samples(errorbar) == samples
    # Test that we recover original errorbar
    assert ev1.get_errorbar(samples) == errorbar
    # Test that half the errorbar needs quadruple samples
    assert ev1.get_samples(errorbar / 2) == 4 * samples
    # Test that twice the samples results in 1/sqrt(2) errorbar
    assert match_to_tol(ev1.get_errorbar(2 * samples), errorbar * 2**-0.5, 0.1)
    # Test that samples is always > 0
    assert ev1.get_samples(1e10) == 1
    # Test the effective variance
    assert match_to_tol(ev1.var_eff, samples * errorbar**2)

    # Test array init
    samples2 = [5, 6]
    errorbar2 = [6.0, 4.0]
    ref = array(samples2) * array(errorbar2)**2
    ev2 = EffectiveVariance(samples2, errorbar2)
    assert ev2.var_eff > min(ref)
    assert ev2.var_eff < max(ref)

    # Test addition of errorbar data
    ev1.add_errorbar_data(6, 4.0)
    assert match_to_tol(ev1.var_eff, ev2.var_eff)

    # Test nodes
    with raises(ValueError):
        ev1.get_samples(0.1, nodes=0)
    # end with
    with raises(ValueError):
        ev1.get_errorbar(5, nodes=0)
    # end with
    nodes = 5
    samples = 15
    errorbar = 0.4
    s1 = ev2.get_samples(errorbar, nodes=1)
    sn = ev2.get_samples(errorbar, nodes=nodes)
    # nodes * samples must approximately match
    assert abs(s1 - nodes * sn) < nodes
    e1 = ev2.get_errorbar(errorbar, nodes=1)
    en = ev2.get_errorbar(errorbar, nodes=nodes)
    # errorbar / sqrt(nodes) must remain constant
    assert match_to_tol(e1, nodes**-0.5 * en)

    # Test addition
    ev3 = EffectiveVariance(7, 3.0)
    ev_add = ev2 + ev3
    assert match_to_tol(ev_add.error_data[:-1], ev2.error_data)
    assert match_to_tol(ev_add.error_data[-1], ev3.error_data[-1])
    assert not match_to_tol(ev_add.var_eff, ev2.var_eff)
    assert not match_to_tol(ev_add.var_eff, ev3.var_eff)

    # Test var_eff init
    var_eff = 4.0
    ev4 = EffectiveVariance(var_eff=var_eff)
    assert match_to_tol(ev4.var_eff, var_eff)
    assert match_to_tol(ev4.get_errorbar(samples=1), var_eff**0.5)
    assert match_to_tol(ev4.get_samples(error=var_eff / 16), 64)
    # Test mean
    ev4.add_var_eff(6.0)
    assert match_to_tol(ev4.var_eff, 5.0)

# end def
