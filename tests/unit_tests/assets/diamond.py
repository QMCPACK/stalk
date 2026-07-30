#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import array

from .helper import morse


# test diamond
a0 = 1.768
h0 = 0.5
pos_diamond = array([
    [0.0, 0.0, 0.0],
    [a0 / 2, a0 / 2, a0 / 2]
])
elem_diamond = ['C', 'C']
axes_diamond = array([
    [a0, a0, 0],
    [0, a0, a0],
    [a0, 0, a0]
])
hessian_diamond = array([[h0]])


def forward_diamond(pos, axes):
    a = axes[0, 0]
    return array([a])
# end def


def backward_H2O(params):
    a = params[0]
    axes = [
        [a, a, 0],
        [0, a, a],
        [a, 0, a]
    ]
    pos = [
        [0.0, 0.0, 0.0],
        [a / 2, a / 2, a / 2]
    ]
    return pos, axes
# end def


def pes_diamond(a, sigma=0.0, path=None):
    V = morse([a0, h0, 0.5, 0.0], a)
    return V, sigma
# end def
