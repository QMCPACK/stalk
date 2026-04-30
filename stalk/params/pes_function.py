#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import isscalar
from scipy.optimize import minimize

from stalk.params.effective_variance import EffectiveVariance
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_result import PesResult
from stalk.util.function_caller import FunctionCaller


class PesFunction(FunctionCaller):

    def evaluate(
        self,
        structure: ParameterSet,
        sigma=0.0,
        add_sigma=False,
        var_eff_map=None,
        interactive=False,  # catch interactive
        **kwargs  # path, dep_jobs
    ):
        result = self._evaluate_structure(
            structure,
            sigma=sigma,
            var_eff_map=var_eff_map,
            **kwargs
        )
        # Load hook to be used in derived classes
        self._load_structure(
            structure,
            result=result,
            sigma=sigma,
            add_sigma=add_sigma,
            var_eff_map=var_eff_map,
        )
    # end def

    def evaluate_all(
        self,
        structures: list[ParameterSet],
        sigmas=None,
        add_sigma=False,
        interactive=False,  # catch interactive
        **kwargs  # path, dep_jobs, var_eff_map
    ):
        if sigmas is None:
            sigmas = len(structures) * [0.0]
        # end if
        for structure, sigma in zip(structures, sigmas):
            self.evaluate(structure, sigma=sigma, add_sigma=add_sigma, **kwargs)
        # end for
    # end def

    def _evaluate_structure(
        self,
        structure: ParameterSet,
        sigma=0.0,
        var_eff_map=None,
        **kwargs
    ):
        eval_args = self.args.copy()
        # Override with kwargs
        eval_args.update(**kwargs)
        # Set samples hook
        self._set_samples(structure, sigma=sigma, var_eff_map=var_eff_map)
        value, error = self.func(structure, sigma=sigma, **eval_args)
        return PesResult(value, error)
    # end def

    def _set_samples(
        self,
        structure: ParameterSet,
        sigma=0.0,
        var_eff_map: EffectiveVarianceMap = None
    ):
        if var_eff_map is not None and sigma > 0.0:
            structure.samples = var_eff_map.get_samples(structure, sigma)
        else:
            structure.samples = None
        # end if
    # end def

    def _update_var_eff_map(
        self,
        structure: ParameterSet,
        var_eff_map: EffectiveVarianceMap
    ):
        if isinstance(var_eff_map, EffectiveVarianceMap) and (
                hasattr(structure, 'samples') and isscalar(structure.samples)):
            # Add the effective variance to the map
            var_eff = EffectiveVariance(structure.samples, structure.error)
            var_eff_map.add_var_eff(structure, var_eff)
        # end if
    # end def

    def _load_structure(
        self,
        structure: ParameterSet,
        result: PesResult = None,
        add_sigma=False,
        sigma=0.0,
        var_eff_map: EffectiveVarianceMap = None
    ):
        if result is None:
            return
        # end if
        # Loading is trivial because the result is readily provided
        if add_sigma:
            result.add_sigma(sigma)
        # end if
        structure.value = result.value
        structure.error = result.error
        self._update_var_eff_map(structure, var_eff_map=var_eff_map)
    # end def

    def relax(
        self,
        structure: ParameterSet,
        **kwargs
    ):
        # Relax numerically using a wrapper around SciPy minimize
        def relax_aux(p):
            s = structure.copy(params=p)
            self.evaluate(s)
            return s.value
        # end def
        p0 = structure.params
        res = minimize(relax_aux, p0, **kwargs)
        structure.set_params(res.x)
    # end def

# end class
