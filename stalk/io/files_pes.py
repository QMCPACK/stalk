#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.io.pes_loader import PesLoader
from stalk.io.xyz_geometry import XyzGeometry
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_function import NotEvaluatedException, PesFunction


def write_xyz_sigma(
    structure: ParameterSet,
    suffix='structure.xyz',
    **kwargs
):
    g = XyzGeometry(suffix=suffix)
    g.write(structure=structure, path=structure.path)
    # Note: creating sigma.dat has been moved to default implementation of the function.
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
        super().__init__(func=func, args=args, **kwargs)
        self.loader = loader
        # Always true for the files PES
        self.create_files = True
    # end def

    def _generate_structure(
        self,
        structure: ParameterSet,
        path='',
        sigma=0.0,
        samples=None,
        var_eff_map: EffectiveVarianceMap = None,
        dep_jobs=None,  # catch dep_jobs
        interactive=False,  # catch interactive
        **kwargs
    ):
        # Store the file path to the structure
        structure.path = self._get_path(structure, path)
        # Associate the sigma with the structure
        structure.sigma = sigma
        # Set the number of samples
        self._set_samples(structure, var_eff_map=var_eff_map, samples=samples)
        # Use params.dat to determine if the jobs have been already generated
        if self.params_file.exists(structure.path):
            print(f'Input files in {structure.path} are already generated. Not regenerating.')
        else:
            self._create_files(structure)
            # Create the jobs and store them in the structure
            # Hot update of eval_args
            eval_args = self.args.copy()
            eval_args.update(**kwargs)
            # Call for the evaluation function
            self.func(structure, sigma=sigma, **eval_args)
        # end if
        # Try to load the value from disk if it exists
        self._try_load_value(structure)
    # end def

    def _evaluate_structure(
        self,
        structure: ParameterSet,
        interactive: bool = False,
        dep_jobs=[],
    ) -> None:
        # Nothing to be done here
        pass
    # end def

    def _finalize_structure(
        self,
        structure: ParameterSet,
        add_sigma: bool = False,
        var_eff_map: EffectiveVarianceMap = None,
        warn_limit=2.0,
        interactive: bool = False,
    ) -> None:
        if structure.valid:
            print(f'Structure {structure.label} is already valid. Not re-evaluating.')
            return
        # end if
        # Try to load the result from disk
        try:
            result = self.loader.load(structure.path)
            if add_sigma:
                result.add_sigma(structure.sigma)
            # end if
            structure.value = result.value
            structure.error = result.error
            # TODO: interactively discard bad data?
            self._warn_energy(structure, warn_limit=warn_limit)
            # Nothing to do here but update the var_eff_map if needed
            self._update_var_eff_map(structure, var_eff_map=var_eff_map)
        except NotEvaluatedException:
            msg = f'Structure {structure.path}/{structure.label} has not been evaluated. '
            msg += 'Supply output file to disk to continue.'
            print(msg)
        # end try
    # end def

# end class
