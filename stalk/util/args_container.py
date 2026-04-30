#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


class ArgsContainer():
    _args = {}

    @property
    def args(self) -> dict:
        return self._args
    # end def

    @args.setter
    def args(self, args):
        if isinstance(args, dict):
            self._args = args
        elif args is None:
            self._args = {}
        else:
            raise TypeError("The argument list must be a dictionary")
        # end if
    # end

    def __init__(self, **args):
        self.args = args
    # end def

    def get_updated(self, new_args: dict) -> dict:
        args = self.args.copy()
        args.update(new_args)
        return args
    # end def

# end class
