#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.io.txt_data import TxtData
from stalk.params.geometry_result import GeometryResult
from stalk.params.parameter_set import ParameterSet
from stalk.params.parameter_structure import ParameterStructure
from stalk.params.pes_function import NotEvaluatedException, PesFunction


class RelaxFunction(PesFunction):
    params_relax_file: TxtData = None

    def __init__(
        self,
        func,
        args: dict = {},  # Keep 'args' for backward compatibility
        **kwargs,  # disable_failed=False, create_files=True, custom kwargs
    ):
        # Init the Pes function
        super().__init__(func, args=args, **kwargs)
        # Initialize the init/relaxed params file handlers
        self.params_file = TxtData('params_init.dat')
        self.params_relax_file = TxtData('params_relax.dat')
    # end def

    # Override the evaluate method to handle relaxation
    def _evaluate_structure(
        self,
        structure: ParameterSet,
        interactive: bool = False,
        reset_value=False,
        dep_jobs=[],
    ) -> None:
        if interactive:
            self._prompt([structure])
        # end if
        if reset_value:
            structure.reset_value()
        # end if
        # Try to load relaxed parameters from disk
        params_relax = self.params_relax_file.load_result(structure.path, None, ndmin=1)
        if params_relax is not None:
            print(f'{structure.path} is already relaxed.')
            value, error = structure.value, structure.error
            structure.params = params_relax
            structure.value, structure.error = value, error
            return
        # end if
        try:
            # Unlike elsewhere, the structure, value and error can be set here directly and
            # bypass using loaders. The data is transferred via the structure
            res = self.func(structure, **self.args)
            if isinstance(res, ParameterSet):
                structure.params = res.params
                structure.value = res.value
                structure.error = res.error
            elif res is None and isinstance(structure, ParameterStructure):
                # Must be a GeometryLoader, providing a GeometryResult
                geom_result = self.loader.load(structure.path)
                if not isinstance(geom_result, GeometryResult):
                    raise TypeError("The loader must return a GeometryResult!")
                # end if
                structure.pos = geom_result.pos
                structure.axes = geom_result.axes
            else:
                raise ValueError("RelaxFunction must return a ParameterSet containing the relaxed parameters, or it must be supplied with a GeometryLoader")
            # end if
            # Ideally, the relax function also evaluates the energy but it is not strictly required.
            if structure.evaluated:
                super()._save_value(structure)
            # end if
            if self.create_files and structure.path is not None:
                self.params_relax_file.save_result(structure.path, structure.params)
            # end if
        except NotEvaluatedException:
            print(f'{structure.path} could not be relaxed.')
        # end try
    # end def

    def _finalize_structure(
        self,
        structure: ParameterSet,
        **kwargs  # add_sigma=False, var_eff_map=None, warn_limit=2.0, interactive=False
    ) -> None:
        # Nothing to do here
        pass
    # end def

    def __call__(
        self,
        structure: ParameterSet,
        **kwargs
    ) -> ParameterSet:
        self.evaluate(structure, **kwargs)
        return structure
    # end def

# end class
