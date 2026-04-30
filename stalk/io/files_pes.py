#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from os import makedirs
import sys

from numpy import isnan
from stalk.io.pes_loader import PesLoader
from stalk.params.parameter_set import ParameterSet
from stalk.params.pes_function import PesFunction
from stalk.util.util import directorize
from stalk.io.util import write_xyz_sigma


class FilesPes(PesFunction):
    loader = None

    def __init__(
        self,
        func=write_xyz_sigma,
        args={},
        loader=None,
        load_args={},
        suffix='energy.dat'
    ):
        self.func = func
        self.args = args
        if isinstance(loader, PesLoader):
            self.loader = loader
        else:
            self.loader = PesLoader(suffix=suffix, **load_args)
        # end if
    # end def

    def evaluate(
        self,
        structure: ParameterSet,
        sigma=0.0,
        add_sigma=False,
        **kwargs  # path, interactive, dep_jobs
    ):
        result = self._evaluate_structure(structure, sigma=sigma, **kwargs)
        # Load hook to be used in derived classes
        finished = self._load_structure(
            structure,
            result=result,
            sigma=sigma,
            add_sigma=add_sigma
        )
        if not finished:
            sys.exit("The job has been finished yet.")
        # end if
    # end def

    def evaluate_all(
        self,
        structures: list[ParameterSet],
        sigmas=None,
        add_sigma=False,
        **kwargs  # path, interactive, dep_jobs
    ):
        if sigmas is None:
            sigmas = len(structures) * [0.0]
        # end if
        finished = True
        for structure, sigma in zip(structures, sigmas):
            result = self._evaluate_structure(structure, sigma=sigma, **kwargs)
            # Load hook to be used in derived classes
            finished &= self._load_structure(
                structure,
                result=result,
                sigma=sigma,
                add_sigma=add_sigma
            )
        # end for
        if not finished:
            sys.exit("Some jobs have not been finished yet.")
        # end if
    # end def

    def _evaluate_structure(
        self,
        structure: ParameterSet,
        path='',
        sigma=0.0,
        dep_jobs=None,  # catch dep_jobs
        interactive=False,  # catch interactive
        **kwargs
    ):
        file_path = f'{directorize(path)}{structure.label}/'
        eval_args = self.args.copy()
        # Override with kwargs
        eval_args.update(**kwargs)
        structure.file_path = file_path
        makedirs(file_path, exist_ok=True)
        self.func(
            structure,
            sigma=sigma,
            **eval_args
        )
    # end def

    def _load_structure(
        self,
        structure: ParameterSet,
        add_sigma=False,
        sigma=0.0,
        # warn_limit=2.0,
        **kwargs
    ):
        result = self.loader.load(structure.file_path)
        if isnan(result.value):
            finished = False
        else:
            if add_sigma and hasattr(structure, "sigma"):
                result.add_sigma(sigma)
            # end if
            structure.value = result.value
            structure.error = result.error
            finished = True
        # end if
        return finished
    # ne dedf

# end class
