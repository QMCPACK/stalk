#!/usr/bin/env python3

from stalk.io.txt_data import TxtData
from stalk.params.pes_function import NotEvaluatedException
from stalk.params.pes_result import PesResult
from stalk.util.args_container import ArgsContainer

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class PesLoader(ArgsContainer, TxtData):
    _suffix = 'energy.dat'

    def __init__(
        self,
        args: dict = {},  # Keep 'args' for backward compatibility
        scale=1.0,
        **kwargs,
    ):
        args.update(**kwargs)
        # suffix = None means we'll use class-level default
        suffix = args.pop('suffix', None)
        TxtData.__init__(self, suffix=suffix, scale=scale)
        ArgsContainer.__init__(self, **args)
    # end def

    def load(self, path: str) -> PesResult:
        # Loading hook
        res = self._load(path, **self.args)
        print(f'Loaded energy from {path}: {res.value} ± {res.error}')
        return res
    # end def

    # Loading hook that can be overridden in derived classes for custom loading behavior
    def _load(self, path: str, **kwargs) -> PesResult:
        try:
            # PesLoader will not tolerate missing files unless default=nan is in kwargs
            data = self.load_result(path, **kwargs)
            if len(data) == 1:
                result = PesResult(data[0])
            else:
                result = PesResult(data[0], data[1])
            # end if
        except FileNotFoundError as e:
            msg = f'PesLoader could not find energy file in {path}. '
            msg += 'To continue, create one with proper energy data or NaN.'
            raise NotEvaluatedException(msg) from e
        except TypeError as e:
            msg = f'PesLoader failed to load the energy output in {path}.'
            raise TypeError(msg) from e
        # end try
        return result
    # end def

# end class
