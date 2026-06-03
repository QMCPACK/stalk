#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from os import makedirs

from stalk.io.pes_loader import PesLoader
from stalk.io.txt_data import TxtData
from stalk.io.xyz_geometry import XyzGeometry
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_function import NotEvaluatedException, PesFunction
from stalk.util.util import directorize


def write_xyz_sigma(
    structure: ParameterSet,
    suffix='structure.xyz',
    sigma=None,
    sigma_suffix='sigma.dat',
    **kwargs
):
    g = XyzGeometry(suffix=suffix)
    g.write(structure=structure, path=structure.path)
    if sigma is not None:
        s = TxtData(suffix=sigma_suffix)
        s.save_result(structure.path, [sigma])
    # end if
# end def


class FilesPes(PesFunction):
    loader: PesLoader = None

    def __init__(
        self,
        func=write_xyz_sigma,
        args={},
        loader: PesLoader = PesLoader(),
        **kwargs  # extra kwargs for PES evaluation
    ):
        # Init the function caller
        super().__init__(func=func, args=args, **kwargs)
        self.loader = loader
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
        # Write the file path to the structure
        structure.path = f'{directorize(path)}{structure.label}/'
        # Hot update of eval_args
        eval_args = self.args.copy()
        eval_args.update(**kwargs)
        # Associate the sigma with the structure
        structure.sigma = sigma
        # Set the number of samples
        self._set_samples(structure, var_eff_map=var_eff_map, samples=samples)
        # Write the input files to disk
        makedirs(structure.path, exist_ok=True)
        # Call for the evaluation function
        self.func(
            structure,
            sigma=sigma,
            **eval_args
        )
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
            print(f'Structure {structure.label} has not be evaluated. Supply output file to disk to continue.')
        # end try
    # end def

# end class
