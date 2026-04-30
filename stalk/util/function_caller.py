#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


from stalk.util.args_container import ArgsContainer


class FunctionCaller(ArgsContainer):
    _func = None

    @property
    def func(self):
        return self._func
    # end def

    @func.setter
    def func(self, func):
        if callable(func):
            self._func = func
        else:
            raise TypeError("The function must be callable!")
        # end if
    # end

    def __init__(
        self,
        func,
        args: dict = {},  # keep positional 'args' in place for backward compatibility
        **kwargs
    ):
        if isinstance(func, FunctionCaller):
            self.func = func.func
            args.update(kwargs)
            self.args = func.get_updated(args)
        else:
            self.func = func
            super().__init__(**args, **kwargs)
        # end if
    # end def

# end class
