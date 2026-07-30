#!/usr/bin/env python

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array
from stalk.util import match_to_tol


def test_parameter_tools():
    from stalk.params import distance, bond_angle, mean_distances
    # water molecule
    pos = array('''
    0.00000        0.00000        0.11779
    0.00000        0.75545       -0.47116
    0.00000       -0.75545       -0.47116
    '''.split(), dtype=float).reshape(-1, 3)
    assert match_to_tol(distance(pos[0], pos[1]), 0.957897074324)
    assert match_to_tol(distance(pos[0], pos[2]), 0.957897074324)
    assert match_to_tol(distance(pos[1], pos[2]), 1.5109)
    assert match_to_tol(mean_distances(
        [(pos[0], pos[1]), (pos[0], pos[2])]), 0.957897074324)
    assert match_to_tol(bond_angle(pos[1], pos[0], pos[2]), 104.1199307245)
# end def
