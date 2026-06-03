#!/usr/bin/env python3
'''A wrapper class for generating Nexus functions to produce and represent a PES.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import isscalar
from pickle import load

from nexus import run_project, bundle

from stalk.io.pes_loader import PesLoader
from stalk.nexus.nexus_structure import NexusStructure
from stalk.params.pes_function import NotEvaluatedException, PesFunction
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.util.util import directorize


class NexusPes(PesFunction):
    loader: PesLoader = None
    disable_failed = False
    bundle_jobs = False

    def __init__(
        self,
        func,
        args: dict = {},  # Keep 'args' for backward compatibility
        loader: PesLoader = None,
        disable_failed=False,
        bundle_jobs=False,
        **kwargs,
    ):
        # Init the function caller
        super().__init__(func, args=args, **kwargs)
        self.disable_failed = disable_failed
        self.bundle_jobs = bundle_jobs
        self.loader = loader
    # end def

    # Override generation function to support Nexus job generation
    def _generate_structure(
        self,
        structure: NexusStructure,
        path='',
        sigma=0.0,
        samples=None,
        var_eff_map: EffectiveVarianceMap = None,
        interactive=False,
        dep_jobs=[],
        # Track repeated Nexus identifiers to avoid redundant job generation
        gen_paths: set = set(),
        **kwargs
    ) -> None:
        # Do not redo jobs
        if structure.generated:
            return
        # end if
        # Write the file path to the structure
        structure.path = f'{directorize(path)}{structure.label}/'
        # Hot update of eval_args
        eval_args = self.args.copy()
        eval_args.update(**kwargs)
        # Associate the sigma with the structure
        structure.sigma = sigma
        # Set the number of samples
        self._set_samples(structure, var_eff_map=var_eff_map, samples=samples)
        # Create the jobs and store them in the structure
        if structure.path in gen_paths:
            structure.jobs = []
        else:
            structure.jobs = self.func(structure, dep_jobs=dep_jobs, **eval_args)
            gen_paths.add(structure.path)
        # end if
    # end def

    def _generate_structure_all(
        self,
        structures: list[NexusStructure],
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
        # Generate the jobs only in unique paths
        gen_paths = set()
        for structure, sigma in zip(structures, sigmas):
            self._generate_structure(
                structure,
                path=path,
                sigma=sigma,
                var_eff_map=var_eff_map,
                interactive=interactive,
                dep_jobs=dep_jobs,
                # track repeated Nexus identifiers
                gen_paths=gen_paths,
                **kwargs
            )
        # end for
    # end def

    def _evaluate_structure(
        self,
        structure: NexusStructure,
        interactive: bool = False,
        dep_jobs=[],
    ) -> None:
        if interactive:
            self._prompt([structure])
        # end if
        # TODO: Try to load result from disk first
        # Run Nexus jobs
        jobs = dep_jobs + structure.jobs
        if self.bundle_jobs:
            run_project(bundle(jobs))
        else:
            run_project(jobs)
        # end if
    # end def

    def _evaluate_structure_all(
        self,
        structures: list[NexusStructure],
        interactive: bool = False,
        dep_jobs=[],
    ) -> None:
        if interactive:
            self._prompt(structures)
        # end if
        jobs = dep_jobs
        for structure in structures:
            jobs += structure.jobs
        # end for
        if self.bundle_jobs:
            run_project(bundle(jobs))
        else:
            run_project(jobs)
        # end if
    # end def

    def _finalize_structure(
        self,
        structure: NexusStructure,
        add_sigma: bool = False,
        var_eff_map: EffectiveVarianceMap = None,
        warn_limit=2.0,
        interactive: bool = False,
    ):
        # Then, try to load the result
        try:
            result = self.loader.load(structure.path)
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
        except NotEvaluatedException as e:
            if self.disable_failed:
                structure.enabled = False
                print(f'Failed to load result for {structure.label} from {structure.path}. Disabling structure.')
            else:
                raise NotEvaluatedException(f'Failed to load result for {structure.label} from {structure.path}.') from e
            # end if
        # end try
    # end def

    def _prompt(self, structures: list[NexusStructure]):
        new_job_strs = []
        for structure in structures:
            if structure.generated:
                for job in structure.jobs:
                    sim_path = '{}/sim_{}/sim.p'.format(job.path, job.identifier)
                    finished = False
                    try:
                        with open(sim_path, mode='rb') as f:
                            sim = load(f)
                            finished = sim.finished
                        # end with
                    except (FileNotFoundError, AttributeError):
                        pass
                    # end try
                    if not finished:
                        job_str = '  {}'.format(job.path)
                        if hasattr(job, "samples") and isscalar(job.samples):
                            job_str += f' ({job.samples}x samples)'
                        # end if
                        new_job_strs.append(job_str)
                    # end if
                # end for
            # end if
        # end for
        if len(new_job_strs) > 0:
            print("About to submit the following jobs:")
            for job_str in new_job_strs:
                print(job_str)
            # end for
            proceed = input("Proceed (Y/n)? ")
            if proceed == 'n':
                exit("Submission cancelled by user.")
            # end if
        # end if
    # end def

    def relax(
        self,
        *args,
        **kwargs
    ):
        msg = "Relaxation not implemented in NexusPes class, use NexusGeometry instead"
        raise NotImplementedError(msg)
    # end def

# end class
