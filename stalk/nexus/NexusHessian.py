#!/usr/bin/env python3
"""NexusHessian class inherits ParameterHessian with Nexus functions.
"""

from numpy import array

from stalk.util import directorize
from stalk.io.PesLoader import PesLoader
from stalk.params.ParameterHessian import ParameterHessian
from stalk.nexus.NexusGenerator import NexusGenerator

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class NexusHessian(ParameterHessian):
    """NexusHessian class inherits ParameterHessian with Nexus functions.
    """

    def _evaluate_energies(
        self,
        structure_list,
        label_list=None,
        path='fdiff',
        pes=None,
        pes_func=None,
        pes_args={},
        loader=None,
        load_func=None,
        load_args={},
        **kwargs,
    ):
        # Generate jobs
        if not isinstance(pes, NexusGenerator):
            # Checks are made in the wrapper class
            pes = NexusGenerator(pes_func, pes_args)
        # end if
        jobs = []
        for s, label in zip(structure_list, label_list):
            dir = '{}{}'.format(directorize(path), label)
            # Make a copy structure for job generation
            jobs += pes.generate(s.copy(), dir)
        # end for
        from nexus import run_project
        run_project(jobs)

        # Load jobs
        if not isinstance(loader, PesLoader):
            # Try to instantiate PesLoader class in the hopes that load_func conforms
            loader = PesLoader(load_args)
            loader.__load__ = load_func
        # end if
        Es = []
        for label in label_list:
            dir = '{}{}'.format(directorize(path), label)
            E = loader.load(path=dir).get_value()
            Es.append(E)
        # end for
        return array(Es)
    # end def

# end class
