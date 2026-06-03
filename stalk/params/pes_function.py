#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings

from numpy import isscalar
from scipy.optimize import minimize

from stalk.params.effective_variance import EffectiveVariance
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_result import PesResult
from stalk.util.function_caller import FunctionCaller
from stalk.util.util import directorize


class PesFunction(FunctionCaller):

    def evaluate(
        self,
        structure: ParameterSet,
        path='',
        sigma=0.0,
        samples=None,
        add_sigma=False,
        var_eff_map: EffectiveVarianceMap = None,
        interactive=False,
        warn_limit=2.0,
        dep_jobs=[],
        **kwargs  # samples, etc.
    ) -> None:
        # Generate hook creates the structure file and sigma
        self._generate_structure(
            structure,
            path=path,
            sigma=sigma,
            samples=samples,
            var_eff_map=var_eff_map,
            interactive=interactive,
            dep_jobs=dep_jobs,
            **kwargs
        )
        # Evaluate hook where the structure is evaluated and the result is set to the structure
        self._evaluate_structure(
            structure,
            interactive=interactive,
            dep_jobs=dep_jobs,
        )
        # Finalization hook where the values are loaded and can be used to update the process
        self._finalize_structure(
            structure,
            add_sigma=add_sigma,
            var_eff_map=var_eff_map,
            warn_limit=warn_limit,
            interactive=interactive,
        )
    # end def

    def evaluate_all(
        self,
        structures: list[ParameterSet],
        sigmas=None,
        path='',
        add_sigma=False,
        var_eff_map=None,
        interactive=False,
        warn_limit=2.0,
        dep_jobs=[],
        **kwargs
    ) -> None:
        # Generate hook creates the structure file and sigma
        self._generate_structure_all(
            structures,
            path=path,
            sigmas=sigmas,
            var_eff_map=var_eff_map,
            dep_jobs=dep_jobs,
            **kwargs
        )
        # Evaluate hook where the structure is evaluated and the result is set to the structure
        self._evaluate_structure_all(
            structures,
            interactive=interactive,
            dep_jobs=dep_jobs,
        )
        # Finalization hook where the values are loaded and can be used to update the process
        self._finalize_structure_all(
            structures,
            add_sigma=add_sigma,
            var_eff_map=var_eff_map,
            warn_limit=warn_limit,
            interactive=interactive,
        )
    # end def

    def _generate_structure(
        self,
        structure: ParameterSet,
        path: str,
        sigma=0.0,
        samples=None,
        var_eff_map: EffectiveVarianceMap = None,
        interactive=False,
        dep_jobs=[],
        **kwargs
    ) -> None:
        # Write the file path to the structure
        structure.path = f'{directorize(path)}{structure.label}/'
        # TODO: write the appropriate files to disk if not already there
        # Associate the sigma with the structure
        structure.sigma = sigma
        self._set_samples(structure, var_eff_map=var_eff_map, samples=samples)
    # end def

    def _generate_structure_all(
        self,
        structures: ParameterSet,
        path,
        sigmas=None,
        var_eff_map=None,
        interactive=False,
        dep_jobs=[],
        **kwargs
    ) -> None:
        if sigmas is None:
            sigmas = [0.0] * len(structures)
        # end if
        # In the default implementation, just call the single structure version for each structure
        for structure, sigma in zip(structures, sigmas):
            self._generate_structure(
                structure,
                path=path,
                sigma=sigma,
                var_eff_map=var_eff_map,
                interactive=interactive,
                dep_jobs=dep_jobs,
                **kwargs
            )
        # end for
    # end def

    def _evaluate_structure(
        self,
        structure: ParameterSet,
        interactive: bool = False,
        dep_jobs=[],
    ) -> None:
        if interactive:
            self._prompt([structure])
        # end if
        # TODO: Try to load result from disk first
        try:
            value, error = self.func(structure, **self.args)
            # Unlike elsewhere, the value and error can be set here directly and bypass
            # using loaders. The data is transferred via the structure
            structure.value = value
            structure.error = error
        except NotEvaluatedException:
            print(f'Structure {structure.label} could not be evaluated.')
        # end try
    # end def

    def _evaluate_structure_all(
        self,
        structures: ParameterSet,
        interactive: bool = False,
        dep_jobs=[],
    ) -> None:
        if interactive:
            self._prompt(structures)
        # end if
        # In the default implementation, just call the single structure version for each structure
        for structure in structures:
            self._evaluate_structure(structure)
        # end for
    # end def

    def _finalize_structure(
        self,
        structure: ParameterSet,
        add_sigma: bool = False,
        var_eff_map: EffectiveVarianceMap = None,
        warn_limit=2.0,
        interactive: bool = False,
    ):
        # Moving value+error back and forth to comply with other implementations
        result = PesResult(structure.value, structure.error)
        if add_sigma:
            result.add_sigma(structure.sigma)
        # end if
        structure.value = result.value
        structure.error = result.error
        # TODO: write to disk
        # TODO: interactively discard bad data?
        self._warn_energy(structure, warn_limit=warn_limit)
        # Nothing to do here but update the var_eff_map if needed
        self._update_var_eff_map(structure, var_eff_map=var_eff_map)
    # end def

    def _finalize_structure_all(
        self,
        structures: ParameterSet,
        add_sigma: bool = False,
        var_eff_map: EffectiveVarianceMap = None,
        warn_limit=2.0,
        interactive: bool = False,
    ) -> None:
        # In the default implementation, just call the single structure version for each structure
        for structure in structures:
            self._finalize_structure(
                structure,
                add_sigma=add_sigma,
                var_eff_map=var_eff_map,
                warn_limit=warn_limit,
                interactive=interactive,
            )
        # end for
    # end def

    def _set_samples(
        self,
        structure: ParameterSet,
        var_eff_map: EffectiveVarianceMap = None,
        samples: float | int | None = None,
    ):
        if isinstance(var_eff_map, EffectiveVarianceMap) and structure.sigma > 0.0:
            structure.samples = var_eff_map.get_samples(structure, structure.sigma)
        else:
            structure.samples = samples
        # end if
    # end def

    def _prompt(self, structures: list[ParameterSet]):
        print("About to evaluate the following structures:")
        for structure in structures:
            print(structure.path)
        # end if
        proceed = input("Proceed (Y/n)? ")
        if proceed in ['n', 'N']:
            exit("Submission cancelled by user.")
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
        structure.params = res.x
    # end def

    def _warn_energy(self, structure: ParameterSet, warn_limit=2.0):
        if (structure.sigma > 0.0 and structure.error / structure.sigma > warn_limit):
            msg = f'The error/sigma for {structure.label} is '
            msg += f'{structure.error}/{structure.sigma} '
            msg += f'{structure.error / structure.sigma * 100.0:.2f}%'
            warnings.warn(msg)
        # end if
    # end def

    def get_var_eff(
        self,
        structure: ParameterSet,
        path='path',
        samples: float | int = 10,
        interactive: bool = False,
    ) -> EffectiveVariance:
        self.evaluate(
            structure,
            path=path,
            interactive=interactive,
            samples=samples,
        )
        var_eff = EffectiveVariance(samples, structure.error)
        return var_eff
    # end def

    def get_var_eff_map(
        self,
        structure: ParameterSet,
        path='path',
        samples=10,
        interactive=False,
    ) -> EffectiveVarianceMap:
        var_eff = self.get_var_eff(
            structure,
            path=path,
            interactive=interactive,
            samples=samples,
        )
        var_eff_map = EffectiveVarianceMap(structure, var_eff=var_eff)
        return var_eff_map
    # end def

# end class


# Exception used to indicate that the energy of a structure has not been evaluated yet.
class NotEvaluatedException(Exception):

    def __init__(self, msg):
        super().__init__(self, msg)
    # end def

# end class
