#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises
from stalk.params import ParameterSet
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.effective_variance import EffectiveVariance
from stalk.util.util import match_to_tol


def test_EffectiveVarianceMap():

    with raises(TypeError):
        # Empty init raises error
        EffectiveVarianceMap()
    # end with

    # Test nominal init
    p_init = ParameterSet([1, 2, 3])
    evm0 = EffectiveVarianceMap(p_init)
    assert len(evm0) == 0
    # Cannot get errorbar without data
    with raises(AssertionError):
        evm0.get_samples(p_init, 0.5)
    # end with

    # a simple model for scaling the effective variance
    def get_scaled_vareff(dp_this, d2s_dp2):
        tot_scale = 1.0
        for dp, d2s in zip(dp_this, d2s_dp2):
            tot_scale += d2s * dp**2
        # end for
        return tot_scale
    # end def

    # Add var_eff information
    base_samples = 10
    base_errorbar = 1.0
    base_var_eff = EffectiveVariance(samples=base_samples, error=base_errorbar)
    d2scale_dp2 = [0.2, 0.5, 2.6]
    offsets0 = [-0.2, 0.3, 0.4]
    offsets1 = [-0.5, 0.5]
    offsets2 = [-0.3, -0.1, 0.05, 0.4]
    # Populate data
    for o0 in offsets0:
        dp = [o0, 0, 0]
        p_this = p_init.copy()
        p_this.shift_params(dp)
        var_eff = base_var_eff.var_eff * get_scaled_vareff(dp, d2scale_dp2)
        ev = EffectiveVariance(var_eff=var_eff)
        evm0.add_var_eff(p_this, ev)
    # end for
    for o1 in offsets1:
        dp = [0, o1, 0]
        p_this = p_init.copy()
        p_this.shift_params(dp)
        var_eff = base_var_eff.var_eff * get_scaled_vareff(dp, d2scale_dp2)
        evm0.add_var_eff(p_this, EffectiveVariance(var_eff=var_eff))
    # end for
    for o2 in offsets2:
        dp = [0, 0, o2]
        p_this = p_init.copy()
        p_this.shift_params(dp)
        var_eff = base_var_eff.var_eff * get_scaled_vareff(dp, d2scale_dp2)
        evm0.add_var_eff(p_this, EffectiveVariance(var_eff=var_eff))
    # end for
    # Check data
    for o0 in offsets0:
        dp = [o0, 0, 0]
        p_this = p_init.copy()
        p_this.shift_params(dp)
        var_eff = base_var_eff.var_eff * get_scaled_vareff(dp, d2scale_dp2)
        ev = EffectiveVariance(var_eff=var_eff)
        errorbar = 0.2
        assert match_to_tol(evm0.get_samples(p_this, errorbar), ev.get_samples(errorbar))
    # end for
    for o1 in offsets1:
        dp = [0, o1, 0]
        p_this = p_init.copy()
        p_this.shift_params(dp)
        var_eff = base_var_eff.var_eff * get_scaled_vareff(dp, d2scale_dp2)
        ev = EffectiveVariance(var_eff=var_eff)
        errorbar = 0.1
        assert match_to_tol(evm0.get_samples(p_this, errorbar), ev.get_samples(errorbar))
    # end for
    for o2 in offsets2:
        dp = [0, 0, o2]
        p_this = p_init.copy()
        p_this.shift_params(dp)
        var_eff = base_var_eff.var_eff * get_scaled_vareff(dp, d2scale_dp2)
        ev = EffectiveVariance(var_eff=var_eff)
        errorbar = 0.6
        assert match_to_tol(evm0.get_samples(p_this, errorbar), ev.get_samples(errorbar))
    # end for


# end def
