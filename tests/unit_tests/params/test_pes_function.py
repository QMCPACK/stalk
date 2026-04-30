#!/usr/bin/env python3

from pytest import raises
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_function import PesFunction
from stalk.util.util import match_to_tol

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_PesFunction():

    # Minimal function for testing vs callable
    def pes_func(s: ParameterSet, arg=0, sigma=0.0):
        return sum(s.params), arg
    # end def

    # Test nominal
    args = {"arg": 5.0}
    pf = PesFunction(func=pes_func, args=args)

    # Test evaluation add_sigma=False
    params = [1.0, 2.0, 3.0]
    sigma = 3.5
    s = ParameterSet(params)
    pf.evaluate(s, add_sigma=False, sigma=sigma)
    assert s.value == sum(params)
    assert match_to_tol(s.error, args['arg'])

    pf.evaluate(s, add_sigma=True, sigma=sigma)
    assert s.value != sum(params)
    assert match_to_tol(s.error, (args['arg']**2 + sigma**2)**0.5)

# end def
