#!/usr/bin/env python3

from pytest import raises
from stalk.util.function_caller import FunctionCaller

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_FunctionCaller():

    # Minimal function for testing vs callable
    def func(x, arg=0):
        return x + arg, x + 1
    # end def

    # Test degraded
    with raises(TypeError):
        # Cannot init empty
        FunctionCaller()
    # end with
    with raises(TypeError):
        # func must be callable
        FunctionCaller(func=[])
    # end with
    with raises(TypeError):
        # args must be dict
        FunctionCaller(pes_func=func, pes_args=[])
    # end with

    # Test nominal
    args = {"arg1": 5, "arg2": 15}
    pf = FunctionCaller(func, **args)
    assert pf.args['arg1'] == args['arg1']
    assert pf.args['arg2'] == args['arg2']
    assert pf.func is func

    # Test copy constructor
    pf_copy = FunctionCaller(pf, arg2=6, arg3=7)
    assert pf_copy.args['arg2'] == 6
    assert pf_copy.args['arg3'] == 7
    assert pf_copy.func is func

# end def
