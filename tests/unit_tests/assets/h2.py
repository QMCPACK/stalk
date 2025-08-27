#!/usr/bin/env python3

from numpy import array
from stalk import ParameterStructure
from stalk import ParameterHessian

from .helper import distance, morse

# test H2 molecule
pos_H2 = array('''
0.00000        0.00000        0.7
0.00000        0.00000       -0.7
'''.split(), dtype=float).reshape(-1, 3)
elem_H2 = 'H H'.split()


def forward_H2(pos, **kwargs):
    r = distance(pos[0], pos[1])
    return [r]
# end def


def backward_H2(params, **kwargs):
    H1 = params[0] * array([0.0, 0.0, 0.5])
    H2 = params[0] * array([0.0, 0.0, -0.5])
    return array([H1, H2])
# end def


def forward_H2_alt(pos, factor=2.0, **kwargs):
    r = distance(pos[0], pos[1]) / factor
    return [r]
# end def


def backward_H2_alt(params, factor=2.0, **kwargs):
    H1 = params[0] * array([0.0, 0.0, 0.5 * factor])
    H2 = params[0] * array([0.0, 0.0, -0.5 * factor])
    return array([H1, H2])
# end def


hessian_H2 = array([[1.0]])


def pes_H2(params):  # inaccurate; for testing
    r, = tuple(params)
    V = 0.0
    V += morse([1.4, 1.17, 0.5, 0.0], r)
    return V
# end def


def alt_pes_H2(params):  # inaccurate; for testing
    r, = tuple(params)
    V = 0.0
    V += morse([1.35, 1.17, 0.6, 0.0], r)
    return V
# end def


def get_structure_H2():
    return ParameterStructure(forward=forward_H2, backward=backward_H2, pos=pos_H2, elem=elem_H2)
# end def


def get_hessian_H2():
    return ParameterHessian(hessian=hessian_H2)
# end def
