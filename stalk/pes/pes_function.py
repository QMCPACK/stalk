#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
from numpy import std
import warnings

from numpy import isscalar
from scipy.optimize import minimize

from stalk.pes.pes_loader import PesLoader
from stalk.io.stalk_path import StalkPath
from stalk.params.effective_variance import EffectiveVariance
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.parameter_set import ParameterSet
from stalk.pes.pes_result import PesResult
from stalk.pes.structure_collection import StructureCollection
from stalk.params.util import NotEvaluatedException
from stalk.util.function_caller import FunctionCaller
from stalk.util.noise import Noise


class PesFunction(FunctionCaller):
    disable_failed = False
    # Loader will not be used in the basic implementation
    loader: PesLoader = None  # Optional loader for loading results from disk

    def __init__(
        self,
        func,
        args: dict = {},  # Keep 'args' for backward compatibility
        loader: PesLoader = None,
        disable_failed=False,
        **kwargs,
    ):
        # Init the function caller
        super().__init__(func, args=args, **kwargs)
        self.disable_failed = disable_failed
        self.loader = loader
    # end def

    def generate(
        self,
        structure: ParameterSet,
        path: Path | str | None = None,
        samples=None,
        var_eff_map: EffectiveVarianceMap = None,
        interactive=False,
        dep_jobs=[],
        **kwargs
    ) -> None:
        self._generate_structure(
            structure,
            path=path,
            samples=samples,
            var_eff_map=var_eff_map,
            interactive=interactive,
            dep_jobs=dep_jobs,
            **kwargs
        )
    # end def

    def evaluate(
        self,
        structure: ParameterSet,
        path: Path | str | None = None,
        samples=None,
        add_sigma: bool | Noise = False,
        var_eff_map: EffectiveVarianceMap = None,
        interactive=False,
        warn_limit=2.0,
        dep_jobs=[],
        reset_value=False,
        **kwargs  # etc.
    ) -> None:
        # Generate hook creates the structure file and sigma
        self._generate_structure(
            structure,
            path=path,
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
            reset_value=reset_value,
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
        path: Path | str | None = None,
        samples=None,
        add_sigma: bool | Noise = False,
        var_eff_map=None,
        interactive=False,
        warn_limit=2.0,
        dep_jobs=[],
        reset_value=False,
        **kwargs
    ) -> None:
        # Generate hook creates the structure file and sigma
        self._generate_structure_all(
            structures,
            path=path,
            samples=samples,
            var_eff_map=var_eff_map,
            dep_jobs=dep_jobs,
            **kwargs
        )
        # Evaluate hook where the structure is evaluated and the result is set to the structure
        self._evaluate_structure_all(
            structures,
            interactive=interactive,
            dep_jobs=dep_jobs,
            reset_value=reset_value,
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

    def _set_path(
        self,
        structure: ParameterSet,
        path: str | Path | None,
        required: bool = False
    ) -> None:
        structure.path = StalkPath(path) / structure.label
        if structure.path is None and required:
            raise ValueError('The path must be specified for the structure.')
        # end if
    # end def

    def _generate_structure(
        self,
        structure: ParameterSet,
        path: Path | str | None,
        samples=None,
        var_eff_map: EffectiveVarianceMap = None,
        interactive=False,
        dep_jobs=[],
        **kwargs
    ) -> None:
        """Set structure path and samples and generate the structure input files to disk."""
        # Set path for the structure
        self._set_path(structure, path, required=False)
        # Set the number of samples
        self._set_samples(structure, var_eff_map=var_eff_map, samples=samples)
        # Save input (if enabled)
        structure.save_input(overwrite=False)
    # end def

    def _generate_structure_all(
        self,
        structures: list[ParameterSet],
        path: Path | str | None,
        samples=None,
        var_eff_map=None,
        interactive=False,
        dep_jobs=[],
        **kwargs
    ) -> None:
        # In the default implementation, just call the single structure version for each structure
        for structure in structures:
            self._generate_structure(
                structure,
                path=path,
                samples=samples,
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
        reset_value=False,
        dep_jobs=[],
    ) -> None:
        """Evaluate the structure."""
        if reset_value:
            # Reset value
            structure.reset_value()
        else:
            # Try to load from disk
            structure.try_load_value()
            # Return if already evaluated or the cache loading succeeded
            if structure.evaluated:
                print(f'{structure.path} is already evaluated.')
                return
            # end if
        # end if
        if interactive:
            self._prompt([structure])
        # end if
        try:
            # The raw PES function is expected to return either a scalar or a tuple of (value, error)
            raw_result = self.func(structure, **self.args)
            if isinstance(raw_result, tuple) and len(raw_result) == 2:
                value, error = raw_result
            elif isscalar(raw_result):
                value, error = float(raw_result), 0.0
            else:
                raise ValueError("The PES function must return a scalar or a tuple of (value, error).")
            # end if
            # Unlike elsewhere, the value and error can be set here directly and bypass
            # using loaders. The data is transferred via the structure.
            structure.value = value
            structure.error = error
        except NotEvaluatedException:
            print(f'{structure.label} could not be evaluated.')
        # end try
    # end def

    def _evaluate_structure_all(
        self,
        structures: list[ParameterSet],
        interactive: bool = False,
        dep_jobs=[],
        reset_value=False,
    ) -> None:
        if interactive:
            self._prompt(structures)
        # end if
        # In the default implementation, just call the single structure version for each structure
        for structure in structures:
            self._evaluate_structure(structure, reset_value=reset_value)
        # end for
    # end def

    def _finalize_structure(
        self,
        structure: ParameterSet,
        add_sigma: bool | Noise = False,
        var_eff_map: EffectiveVarianceMap = None,
        warn_limit=2.0,
        interactive: bool = False,
    ):
        # Moving value+error back and forth to comply with other implementations
        result = PesResult(structure.value, structure.error)
        result.add_sigma(structure.sigma, kind=add_sigma)
        structure.value = result.value
        structure.error = result.error
        if structure.path is not None:
            print(f'{structure.path} evaluated to {structure.value:.6f} +/- {structure.error:.6f}.')
        # end if
        # TODO: interactively discard bad data?
        self._warn_energy(structure, warn_limit=warn_limit)
        # Nothing to do here but update the var_eff_map if needed
        self._update_var_eff_map(structure, var_eff_map=var_eff_map)
        # Write to disk (if enabled)
        structure.save_value(overwrite=False)
    # end def

    def _finalize_structure_all(
        self,
        structures: list[ParameterSet],
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
        if var_eff_map is None:
            return
        # end if
        if structure.error > 0.0 and not var_eff_map.empirical:
            var_eff = EffectiveVariance(structure.samples, structure.error)
            var_eff_map.add_var_eff(structure, var_eff)
        # end if
    # end def

    def relax(
        self,
        structure: ParameterSet,
        path=None,
        **kwargs
    ) -> ParameterSet:
        self._set_path(structure, path, required=False)
        params_relax = structure.load(path=structure.path, key='params_out')
        if params_relax is not None:
            structure.params = params_relax
            structure.try_load_value()
            print(f'Loaded relaxed parameters from {structure.path}.')
            return structure
        # end if

        # Relax numerically using a wrapper around SciPy minimize
        def relax_aux(p):
            s = structure.copy(params=p)
            self.evaluate(s, reset_value=True)
            return s.value
        # end def
        p0 = structure.params
        res = minimize(relax_aux, p0, **kwargs)
        structure.params = res.x
        structure.value = relax_aux(res.x)
        # Save relaxed parameters and value to disk if enabled
        structure.save(path=structure.path, overwrite=False, params_out=structure.params)
        structure.save_value(overwrite=False)
        return structure
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
        path=None,
        samples: float | int = 10,
        interactive: bool = False,
        resamples: int | None = None,
    ) -> EffectiveVariance:
        if resamples is None:
            # Use error of the PES evaluate
            self.evaluate(
                structure,
                path=path,
                interactive=interactive,
                samples=samples,
            )
            var_eff = EffectiveVariance(samples, structure.error)
        elif resamples > 1:
            # Use resampling to estimate the effective variance
            structures: list[ParameterSet] = []
            for i in range(resamples):
                s = structure.copy(label=f'resample_{i}')
                structures.append(s)
            # end for
            self.evaluate_all(
                structures,
                path=path,
                interactive=interactive,
                samples=samples,
            )
            error = std([v for v in [s.value for s in structures]])
            var_eff = EffectiveVariance(samples, error, empirical=True)
        else:
            raise ValueError("Resamples must be None or greater than 1.")
        # end if
        return var_eff
    # end def

    def get_var_eff_map(
        self,
        structure: ParameterSet,
        path=None,
        samples=10,
        interactive=False,
        resamples=None,
    ) -> EffectiveVarianceMap:
        var_eff = self.get_var_eff(
            structure,
            path=path,
            interactive=interactive,
            samples=samples,
            resamples=resamples,
        )
        var_eff_map = EffectiveVarianceMap(structure, var_eff=var_eff)
        return var_eff_map
    # end def

    def __call__(
        self,
        structure: ParameterSet | list[ParameterSet] | StructureCollection,
        **kwargs
    ) -> PesResult | list[PesResult]:
        if isinstance(structure, list):
            # Evaluate a list of structures
            self.evaluate_all(structure, **kwargs)
            result = [PesResult(s.value, s.error) for s in structure]
        elif isinstance(structure, ParameterSet):
            # Evaluate a single structure
            self.evaluate(structure, **kwargs)
            result = PesResult(structure.value, structure.error)
        elif isinstance(structure, StructureCollection):
            # Evaluate a list of structures
            result = self.evaluate_all(structure.collect_enabled(), **kwargs)
            # Finalization hook for collective analysis, e.g., line-search
            structure.finalize()
        else:
            raise TypeError("The structure must be a ParameterSet or a list/collection of ParameterSets.")
        # end if
        return result
    # end def

# end class
