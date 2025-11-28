#!/usr/bin/env python

from numpy import array, flipud
from pytest import raises
from stalk.ls.PolynomialFit import PolynomialFit
from stalk.pls.TargetParallelLineSearch import TargetParallelLineSearch
from stalk.util import match_to_tol

from ..assets.h2o import pes_H2O, get_structure_H2O, get_hessian_H2O
from ..assets.helper import Gs_N200_M7

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


# test TargetParallelLineSearch class
def test_TargetParallelLineSearch():

    # test empty init
    with raises(TypeError):
        TargetParallelLineSearch()
    # end with
    structure = get_structure_H2O()
    hessian = get_hessian_H2O()
    srg = TargetParallelLineSearch(
        fit_kind='pf3',
        pes_func=pes_H2O,
        structure=structure,
        hessian=hessian,
        window_frac=0.2,
    )
    assert srg.setup
    assert srg.evaluated
    assert not srg.optimized

    # Test optimization to windows, noises
    windows = [0.11, 0.12]
    noises = [0.013, 0.014]
    M = 5
    N = 10
    srg.optimize(
        fit_kind='pf2',
        windows=windows,
        noises=noises,
        M=M,
        N=N,
    )
    assert srg.optimized
    for tls, W, sigma in zip(srg.ls_list, windows, noises):
        assert tls.optimized
        assert tls.sigma_opt == sigma
        assert tls.W_opt == W
        assert tls.M == M
    # end for
    assert match_to_tol(srg.W_opt, windows)
    assert match_to_tol(srg.sigma_opt, noises)
    assert match_to_tol(srg.M, [M, M])
    assert all(srg.error_d > 0.0)
    assert all(srg.error_p > 0.0)
    assert srg.epsilon_d is None
    assert srg.epsilon_p is None
    assert srg.temperature is None
    statcost_ref = M * sum(array(noises)**-2)
    assert match_to_tol(srg.statistical_cost, statcost_ref)

    # Test optimization to epsilon_d
    epsilon_d = [0.01, 0.02]
    srg.optimize(
        fit_kind='pf3',
        epsilon_d=epsilon_d,
        Gs=[Gs_N200_M7, flipud(Gs_N200_M7)],
        bias_order=1,
        bias_mix=0.0,
    )
    assert srg.optimized
    # Hard-coded references are not externally validated
    windows_ref1 = [0.058189699596450206, 0.13112494623501014]
    noises_ref1 = [0.003879313306430014, 0.005736716397781694]
    assert match_to_tol(srg.W_opt, windows_ref1)
    assert match_to_tol(srg.sigma_opt, noises_ref1)
    for tls, W, sigma in zip(srg.ls_list, windows_ref1, noises_ref1):
        assert tls.optimized
        assert tls.sigma_opt == sigma
        assert tls.W_opt == W
        assert tls.M == 7
        assert tls.target_settings.bias_mix == 0.0
        assert tls.target_settings.bias_order == 1
        # Test if 'pf3' was chosen as the default
        assert tls.target_settings.fit_func == PolynomialFit(3)
    # end for
    assert match_to_tol(srg.M, [7, 7])
    assert match_to_tol(srg.epsilon_d, epsilon_d)
    assert all(srg.error_d < epsilon_d)
    assert all(srg.error_d > 0.0)
    assert all(srg.error_p > 0.0)
    assert srg.epsilon_p is None
    assert srg.temperature is None
    statcost_ref = 7 * sum(array(noises_ref1)**-2)
    assert match_to_tol(srg.statistical_cost, statcost_ref)

    # Test thermal optimization to epsilon_p
    epsilon_p = [0.02, 0.1]
    srg.optimize(
        fit_kind='pf4',
        epsilon_p=epsilon_p,
        # M = 7, N = 100
        Gs=[Gs_N200_M7[:100], flipud(Gs_N200_M7[:100])],
        bias_order=2,
        bias_mix=0.0,
        thermal=True,
    )
    assert srg.optimized
    # Hard-coded references are not externally validated
    windows_ref2 = [0.058189699596450206, 0.13112494623501014]
    noises_ref2 = [0.006465522177383356, 0.009834370967625761]
    assert match_to_tol(srg.W_opt, windows_ref2)
    assert match_to_tol(srg.sigma_opt, noises_ref2)
    for tls, W, sigma in zip(srg.ls_list, windows_ref2, noises_ref2):
        assert tls.optimized
        assert tls.sigma_opt == sigma
        assert tls.W_opt == W
        assert tls.M == 7
        assert tls.target_settings.N == 100
    # end for
    assert match_to_tol(srg.M, [7, 7])
    assert all(srg.error_d > 0.0)
    assert all(srg.error_p > 0.0)
    assert all(srg.error_p < epsilon_p)
    assert all(array(srg.epsilon_d) > 0.0)
    assert match_to_tol(srg.epsilon_p, epsilon_p)
    assert srg.temperature > 0.0
    statcost_ref2 = 7 * sum(array(noises_ref2)**-2)
    assert match_to_tol(srg.statistical_cost, statcost_ref2)

    # Test LS optimization to epsilon_p
    srg = TargetParallelLineSearch(
        fit_kind='pf3',
        pes_func=pes_H2O,
        structure=structure,
        hessian=hessian,
        window_frac=0.2,
    )
    epsilon_p2 = [0.01, 0.03]
    srg.bracket_target_biases()
    srg.optimize(
        fit_kind='pf3',
        epsilon_p=epsilon_p2,
        Gs=[Gs_N200_M7[:150], flipud(Gs_N200_M7[:150])],
        bias_order=2,
        bias_mix=0.0,
        thermal=False,
    )
    assert srg.optimized
    # Hard-coded references are not externally validated
    windows_ref3 = [0.04849141633037517, 0.12292963709532201]
    noises_ref3 = [0.003232761088691678, 0.006556247311750508]
    assert match_to_tol(srg.W_opt, windows_ref3)
    assert match_to_tol(srg.sigma_opt, noises_ref3)
    for tls, W, sigma in zip(srg.ls_list, windows_ref3, noises_ref3):
        assert tls.optimized
        assert tls.sigma_opt == sigma
        assert tls.W_opt == W
        assert tls.M == 7
        assert tls.target_settings.N == 150
    # end for
    assert match_to_tol(srg.M, [7, 7])
    assert all(srg.error_d > 0.0)
    assert all(srg.error_p > 0.0)
    assert all(srg.error_p < epsilon_p)
    assert all(array(srg.epsilon_d) > 0.0)
    assert match_to_tol(srg.epsilon_p, epsilon_p2)
    assert srg.temperature is None
    statcost_ref3 = 7 * sum(array(noises_ref3)**-2)
    assert match_to_tol(srg.statistical_cost, statcost_ref3)

# end def
