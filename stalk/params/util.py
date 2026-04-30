#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import linalg, pi, arccos, array, dot, sin, cos, linspace
from scipy.optimize import minimize

from stalk.params.parameter_set import ParameterSet


def distance(r0, r1):
    '''Return Euclidean distance between two positions'''
    r = linalg.norm(r0 - r1)
    return r
# end def


def nearest_neighbor(r0, r1, axes):
    '''Return nearest neighbor distance and displacement considering periodicity'''
    dr_min = 1e100
    v_min = None
    for v1 in [axes[0], -axes[0], 0]:
        for v2 in [axes[1], -axes[1], 0]:
            for v3 in [axes[2], -axes[2], 0]:
                dr = distance(r0, r1 + v1 + v2 + v3)
                if dr < dr_min:
                    dr_min = dr
                    v_min = v1 + v2 + v3
                # end if
            # end for
        # end for
    # end for
    return dr_min, v_min
# end def


def periodic_distance(r0, r1, axes):
    '''Return minimum distance between two positions considering periodicity'''
    dr_min = nearest_neighbor(r0, r1, axes)
    return dr_min[0]
# end def


def bond_angle(r0, rc, r1, units='ang'):
    '''Return dihedral angle between 3 bodies'''
    v0 = r0 - rc
    v1 = r1 - rc
    ang = angle(v0, v1, units=units)
    return ang
# end def


def periodic_bond_angle(r0, rc, r1, axes, units='ang'):
    '''Return dihedral angle between 3 bodies'''
    r0_nearest = r0 + nearest_neighbor(r0, rc, axes)[1]
    r1_nearest = r1 + nearest_neighbor(r1, rc, axes)[1]
    v0 = r0_nearest - rc
    v1 = r1_nearest - rc
    ang = angle(v0, v1, units=units)
    return ang
# end def


def angle(v0, v1, units='ang'):
    cosang = dot(v0, v1) / linalg.norm(v0) / linalg.norm(v1)
    ang = arccos(cosang) * 180 / pi if units == 'ang' else arccos(cosang)
    return ang
# end def


def mean_distances(pairs, tol=1e-6, axes=None):
    '''Return average distance over (presumably) identical position pairs'''
    rs = []
    for pair in pairs:
        if axes is not None:
            rs.append(periodic_distance(pair[0], pair[1], axes=axes))
        else:
            rs.append(distance(pair[0], pair[1]))
        # end if
    # end for
    return mean_param(rs, tol=tol)
# end def


def mean_param(params, tol=1e-6):
    avg = array([params]).mean()
    if not all(params - avg < tol):
        print("Warning! Some of symmetric parameters stand out:")
        print(params)
    # end if
    return avg
# end def


def invert_pos(pos0, params, forward=None, tol=1.0e-7, method='BFGS'):
    assert forward is not None, 'Must provide forward mapping'

    def dparams_sum(pos1):
        return sum((params - forward(pos1))**2)
    # end def
    pos1 = minimize(dparams_sum, pos0, tol=tol, method=method).x
    return pos1
# end def


def rotate_2d(arr_2d, ang, units='ang'):
    if units == 'ang':
        ang *= pi / 180
    # end if
    R = array([
        [cos(ang), -sin(ang)],
        [sin(ang), cos(ang)]
    ])
    return R @ arr_2d
# end def


def interpolate_params(structure_a: ParameterSet, structure_b: ParameterSet, num_int):
    scales = linspace(0.0, 1.0, num_int + 2)
    dparams = structure_b.params - structure_a.params
    traj = []
    for scale in scales:
        new_params = structure_a.params + scale * dparams
        structure = structure_a.copy(params=new_params)
        traj.append(structure)
    # end for
    return traj
# end def
