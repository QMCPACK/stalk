#!/usr/bin/env python3
'''A wrapper class for generating Nexus functions to produce and represent a PES.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

import warnings
from numpy import isnan, isscalar
from pickle import load

from nexus import run_project, bundle

from stalk.io.pes_loader import PesLoader
from stalk.nexus.nexus_structure import NexusStructure
from stalk.params.pes_function import PesFunction
from stalk.params.pes_result import PesResult
from stalk.params.effective_variance import EffectiveVariance
from stalk.params.effective_variance_map import EffectiveVarianceMap
from stalk.util.util import FF, FP, directorize


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
        super().__init__(func, args=args, **kwargs)
        self.disable_failed = disable_failed
        self.bundle_jobs = bundle_jobs
        self.loader = loader
    # end def

    # Override evaluation function to support job submission and analysis
    def evaluate(
        self,
        structure: NexusStructure,
        sigma=0.0,
        add_sigma=False,
        path='',
        dep_jobs=[],
        interactive=False,
        warn_limit=2.0,
        var_eff_map=None,
        **kwargs
    ):
        # TODO: try to load first, to assess whether to regenerate or not
        self._evaluate_structure(
            structure,
            path=path,
            sigma=sigma,
            dep_jobs=dep_jobs,
            var_eff_map=var_eff_map,
            **kwargs
        )
        if interactive:
            self._prompt([structure])
        # end if
        jobs = dep_jobs + structure.jobs
        run_project(jobs)
        self._load_structure(
            structure,
            add_sigma=add_sigma,
            warn_limit=warn_limit,
            var_eff_map=var_eff_map,
        )
    # end def

    # Override evaluation function to support parallel job submission and analysis
    def evaluate_all(
        self,
        structures: list[NexusStructure],
        sigmas=None,
        add_sigma=False,
        path='',
        interactive=False,
        dep_jobs=[],
        warn_limit=2.0,
        var_eff_map=None,
        **kwargs
    ):
        if sigmas is None:
            sigmas = len(structures) * [0.0]
        # end if
        if dep_jobs is not None:
            jobs = dep_jobs
        else:
            jobs = []
        # end if
        eqm_generated = False
        for structure, sigma in zip(structures, sigmas):
            skip_gen = False
            if structure.label == 'eqm':
                if eqm_generated:
                    skip_gen = True
                else:
                    eqm_generated = True
                    skip_gen = False
                # end if
            # end if
            if not structure.analyzed:
                self._evaluate_structure(
                    structure,
                    path=path,
                    sigma=sigma,
                    dep_jobs=dep_jobs,
                    skip_gen=skip_gen,
                    var_eff_map=var_eff_map,
                    **kwargs,
                )
                if structure.jobs is not None:
                    jobs += structure.jobs
                # end if
            # end if
        # end for
        # TODO: try to load first, to assess whether to regenerate or not
        if interactive:
            self._prompt(structures)
        # end if
        if self.bundle_jobs:
            run_project(bundle(jobs))
        else:
            run_project(jobs)
        # end if

        # Then, load
        for structure in structures:
            self._load_structure(
                structure,
                add_sigma=add_sigma,
                warn_limit=warn_limit,
                var_eff_map=var_eff_map,
            )
        # end for
    # end def

    def _evaluate_structure(
        self,
        structure: NexusStructure,
        path='',
        sigma=0.0,
        skip_gen=False,
        var_eff_map=None,
        **kwargs
    ):
        file_path = f'{directorize(path)}{structure.label}/'
        eval_args = self.args.copy()
        structure.file_path = file_path
        structure.sigma = sigma
        # Do not redo jobs
        if structure.generated:
            return
        # end if
        # Set samples hook
        self._set_samples(structure, sigma=sigma, var_eff_map=var_eff_map)
        # Override with kwargs
        eval_args = self.get_updated(kwargs)
        if not skip_gen:
            jobs = self.func(
                structure.get_nexus_structure(),
                file_path,
                sigma=sigma,
                **eval_args
            )
            structure.jobs = jobs
        # end if
    # end def

    def _load_structure(
        self,
        structure: NexusStructure,
        add_sigma=False,
        warn_limit=2.0,
        var_eff_map=None,
    ):
        if add_sigma:
            result = self.loader.load(structure.file_path, sigma=structure.sigma)
        else:
            result = self.loader.load(structure.file_path)
        # end if
        self._warn_energy(structure, result, warn_limit=warn_limit)
        # Treat failure
        if isnan(result.value) and self.disable_failed:
            structure.enabled = False
        # end if
        structure.value = result.value
        structure.error = result.error
        # Update var_eff_map hook
        self._update_var_eff_map(structure, var_eff_map)
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

    def _warn_energy(self, structure: NexusStructure, result: PesResult, warn_limit=2.0):
        if (structure.sigma is not None and structure.sigma > 0.0 and result.error / structure.sigma > warn_limit):
            msg = f'The error/sigma for {structure.label} is '
            msg += f'{FF.format(result.error)}/{FF.format(structure.sigma)}'
            msg += f'{FP.format(result.error / structure.sigma * 100)}'
            warnings.warn(msg)
        # end if
    # end def

    def get_var_eff(
        self,
        structure: NexusStructure,
        path='path',
        samples=10,
        interactive=False,
    ):
        self.evaluate(
            structure,
            path=path,
            sigma=None,
            samples=samples,
            interactive=interactive,
        )
        var_eff = EffectiveVariance(samples, structure.error)
        return var_eff
    # end def

    def get_var_eff_map(
        self,
        structure: NexusStructure,
        path='path',
        samples=10,
        interactive=False,
    ):
        var_eff = self.get_var_eff(
            structure,
            path=path,
            samples=samples,
            interactive=interactive
        )
        var_eff_map = EffectiveVarianceMap(structure, var_eff=var_eff)
        return var_eff_map

    def relax(
        self,
        *args,
        **kwargs
    ):
        msg = "Relaxation not implemented in NexusPes class, use NexusGeometry instead"
        raise NotImplementedError(msg)
    # end def

# end class
