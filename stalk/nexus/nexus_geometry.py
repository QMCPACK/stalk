#!/usr/bin/env python3
'''A wrapper class for generating Nexus functions to produce and represent a PES.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path
import warnings

from nexus import run_project

from stalk.io.geometry_loader import GeometryLoader
from stalk.nexus.nexus_structure import NexusStructure
from stalk.pes.relax_function import RelaxFunction


class NexusGeometry(RelaxFunction):
    loader: GeometryLoader = None

    def __init__(
        self,
        func,
        args={},  # keep positional 'args' in place for backward compatibility
        create_files=True,  # NexusGeometry must create files
        **kwargs  # loader=GeometryLoader(), disable_failed=False, custom kwargs
    ):
        super().__init__(func, args=args, create_files=True, **kwargs)
    # end def

    def _generate_structure(
        self,
        structure: NexusStructure,
        path: str | Path,
        interactive=False,
        dep_jobs=[],
        **kwargs
    ):
        # Store the file path to the structure
        self._set_path(structure, path, required=True)
        # Use params.dat to determine if the jobs have been already generated
        if self.params_file.exists(structure.path):
            print(f'Nexus jobs in {structure.path} are already generated. Not regenerating.')
            structure.jobs = []  # Set empty jobs to indicate that the structure is generated
        else:
            self._create_files(structure)
            # Create the jobs and store them in the structure
            # Hot update of eval_args
            eval_args = self.args.copy()
            eval_args.update(**kwargs)
            structure.jobs = self.func(structure, dep_jobs=dep_jobs, **eval_args)
        # end if
    # end def

    def _evaluate_structure(
        self,
        structure: NexusStructure,
        interactive: bool = False,
        dep_jobs=[],
        reset_value=True,
    ) -> None:
        if interactive:
            self._prompt([structure])
        # end if
        # Run Nexus jobs
        jobs = dep_jobs + structure.jobs
        run_project(jobs)
    # end def

    def _finalize_structure(
        self,
        structure: NexusStructure,
        add_sigma: bool = False,
        interactive: bool = False,
        **kwargs,
    ):
        # Try to load relaxed parameters from disk
        params_relax = self.params_relax_file.load_result(structure.path, None, ndmin=1)
        if params_relax is not None:
            print(f'{structure.path} is already relaxed.')
            structure.params = params_relax
            return
        # Then, try to load the result
        res = self.loader.load(structure.path)
        if res.get_pos() is not None:
            structure.set_position(res.get_pos(), res.get_axes())
        else:
            warnings.warn("Running or loading of the relaxation result was unsuccessful!", UserWarning)
        # end if
        if self.create_files and structure.path is not None:
            self.params_relax_file.save_result(structure.path, structure.params)
        # end if
    # end def

# end class
