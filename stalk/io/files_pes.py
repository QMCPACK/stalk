#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path

from stalk.io.xyz_geometry import XyzGeometry
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.parameter_set import ParameterSet
from stalk.pes.pes_function import NotEvaluatedException, PesFunction
from stalk.pes.pes_loader import PesLoader
from stalk.util.noise import Noise


def write_xyz_sigma(
    structure: ParameterSet,
    suffix='structure.xyz',
    **kwargs
):
    g = XyzGeometry(suffix=suffix)
    g.write(structure=structure, path=structure.path)
    # Note: creating sigma.in has been moved to default implementation of the function.
# end def


class FilesPes(PesFunction):

    def __init__(
        self,
        func=write_xyz_sigma,
        args={},
        loader: PesLoader = PesLoader(),
        **kwargs  # disable_failed=False, ...
    ):
        # Init the function caller
        super().__init__(func=func, args=args, loader=loader, **kwargs)
    # end def

    def _generate_structure(
        self,
        structure: ParameterSet,
        path: str | Path,
        samples=None,
        var_eff_map: EffectiveVarianceMap = None,
        dep_jobs=None,  # catch dep_jobs
        interactive=False,  # catch interactive
        **kwargs
    ):
        # Store the file path to the structure
        self._set_path(structure, path, required=True)
        # Set the number of samples
        self._set_samples(structure, var_eff_map=var_eff_map, samples=samples)
        # Use params.in to determine if the jobs have been already generated
        params = structure.load(key='params')
        if params is not None:
            print(f'Input files in {structure.path} are already generated. Not regenerating.')
        else:
            # Create the jobs and store them in the structure
            # Hot update of eval_args
            eval_args = self.args.copy()
            eval_args.update(**kwargs)
            # Call for the evaluation function
            self.func(structure, **eval_args)
            # Save input
            structure.save_input(overwrite=False)
        # end if
    # end def

    def _evaluate_structure(
        self,
        structure: ParameterSet,
        interactive: bool = False,
        dep_jobs=[],
        reset_value: bool = False,
    ) -> None:
        # Nothing to be done here, except assume that value/error be written to disk by
        # the user.
        if reset_value:
            structure.reset_value()
        else:
            structure.try_load_value()
        # end if
    # end def

    def _finalize_structure(
        self,
        structure: ParameterSet,
        add_sigma: bool | Noise = False,
        var_eff_map: EffectiveVarianceMap = None,
        warn_limit=2.0,
        interactive: bool = False,
    ) -> None:
        if structure.valid:
            print(f'{structure.path} is already evaluated. Not re-evaluating.')
            return
        # end if
        # Try to load the result from disk
        try:
            result = self.loader.load(structure.path)
            result.add_sigma(structure.sigma, kind=add_sigma)
            structure.value = result.value
            structure.error = result.error
            # TODO: interactively discard bad data?
            self._warn_energy(structure, warn_limit=warn_limit)
            # Nothing to do here but update the var_eff_map if needed
            self._update_var_eff_map(structure, var_eff_map=var_eff_map)
            # Write to disk
            structure.save_value(overwrite=False)
        except NotEvaluatedException:
            msg = f'{structure.path} has not been evaluated. '
            msg += 'Supply output file to disk to continue.'
            print(msg)
        # end try
    # end def

# end class
