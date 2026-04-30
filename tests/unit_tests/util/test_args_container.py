#!/usr/bin/env python3

from pytest import raises
from stalk.util.args_container import ArgsContainer

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_ArgsContainer():

    # Cannot initialize with positional arguments
    with raises(TypeError):
        ArgsContainer(3)
    # end with

    # Test nominal
    args = {"arg1": 5, "arg2": 15}
    pf = ArgsContainer(**args)
    assert pf.args == args

    # Test getting updated args
    updates = {"arg1": 8, "arg3": 4}
    new_args = pf.get_updated(updates)
    assert all([new_args[key] == updates[key] for key in updates.keys()])
    assert new_args['arg2'] == args['arg2']

    # Test resetting to empty dict
    pf.args = None
    assert len(pf.args) == 0

    with raises(TypeError):
        pf.args = 1
    # end with
# end def
