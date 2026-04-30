#!/usr/bin/env python3

import warnings
from numpy import nan

from stalk.io.util import load_energy
from stalk.params.pes_result import PesResult
from stalk.util.args_container import ArgsContainer
from stalk.util.util import check_result_file

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class PesLoader(ArgsContainer):

    def load(
        self,
        path: str,
        sigma=0.0,
        only_warn=True,  # Only warn instead of raising exception
        **kwargs
    ) -> PesResult:
        # Hot update of args
        args = self.get_updated(kwargs)
        scale = args.pop('scale', 1.0)

        try:
            filename = check_result_file(path, args)
        except FileNotFoundError as e:
            if only_warn:
                warnings.warn(e.args[0])
                return PesResult(nan)
            else:
                raise e.add_note('Aborting.')
            # end if
        # end try
        res = self._load(filename, **args)
        # end if
        # Rescale to model units
        res.rescale(scale)
        # If a non-zero, artificial errorbar is requested, add it to result
        res.add_sigma(sigma)
        print(f'Loaded energy from {filename}')
        return res
    # end def

    def _load(self, structure, **kwargs) -> PesResult:
        res = load_energy(structure)
        return res
    # end def

# end class
