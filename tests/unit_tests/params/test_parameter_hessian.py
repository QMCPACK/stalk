#!/usr/bin/env python

from pytest import raises, warns
from numpy import linalg, where

from stalk.params.pes_function import PesFunction
from stalk.params.parameter_hessian import ParameterHessian

from stalk.util.util import match_to_tol
from unit_tests.assets.h2o import hessian_H2O, get_structure_H2O, pes_H2O

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_ParameterHessian():

    # Test empty
    h0 = ParameterHessian()
    assert h0.structure is None
    assert h0.hessian is None
    assert h0.U is None
    assert h0.directions is None
    assert h0.lambdas is None
    assert h0.require_consistent
    assert h0.disable_limit == 0.0

    with raises(AssertionError):
        h0.hessian = [0]
    # end with

    # Test based on structure only
    s1 = get_structure_H2O()
    h1_ref = [[1., 0.], [0., 1.]]
    U1_ref = [[1., 0.], [0., 1.]]
    Lambda1_ref = [1., 1.]
    h1 = ParameterHessian(structure=s1, disable_limit=0.1)
    assert h1.structure is s1
    assert match_to_tol(h1.hessian, h1_ref)
    assert match_to_tol(h1.U, U1_ref)
    assert match_to_tol(h1.directions, U1_ref)
    assert match_to_tol(h1.lambdas, Lambda1_ref)
    assert match_to_tol(h1.enabled_directions, U1_ref)
    assert match_to_tol(h1.enabled_lambdas, Lambda1_ref)

    # Test based on structure and Hessian
    s2 = s1.copy()
    h2_ref = hessian_H2O
    Lambda2_ref, U2_ref = linalg.eig(hessian_H2O)
    h2 = ParameterHessian(structure=s2, hessian=hessian_H2O)
    assert h2.structure is s2
    assert match_to_tol(h2.hessian, h2_ref)
    assert match_to_tol(h2.U, U2_ref)
    assert match_to_tol(h2.directions, U2_ref.T)
    assert match_to_tol(h2.lambdas, Lambda2_ref)

    # Test computation by finite-difference
    h3_ref = [[4.0, 0.0], [0.0, 1.0]]  # see def pes_H2O()
    E3_ref = -0.5  # see def pes_H2O()
    h3 = ParameterHessian()
    pes = PesFunction(pes_H2O)
    assert h3.hessian is None
    h3.compute_fdiff(
        structure=s1.copy(),
        pes=pes,
        dp=[0.01, 0.02]
    )
    assert match_to_tol(h3.hessian, h3_ref, 1e-3)
    assert match_to_tol(h3.structure.value, E3_ref, 1e-6)

    # Test disabled
    LambdaD_ref, UD_ref = linalg.eig(hessian_H2O)
    # This limit will cut out the lower direction
    hD = ParameterHessian(disable_limit=0.3)
    assert hD.disable_limit == 0.3
    hD.structure = s2
    hD.hessian = hessian_H2O
    idx = where(LambdaD_ref / max(LambdaD_ref) > 0.3)
    assert match_to_tol(hD.enabled_directions, UD_ref.T[idx, :])
    assert match_to_tol(hD.enabled_lambdas, LambdaD_ref[idx])

    # Test warning and scalar dp
    s4 = s1.copy()
    s4.shift_params([0.1, 0.1])
    h4 = ParameterHessian()
    with warns(UserWarning):
        # Expect the energy of eqm_p0-0.01_p1-0.01 to be the lowest -> warning
        h4.compute_fdiff(
            structure=s4,
            pes=pes,
            dp=0.01
        )
    # end with

# end def
