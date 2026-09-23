#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.pes.geometry_result import GeometryResult
from stalk.params.parameter_set import ParameterSet
from stalk.params.parameter_structure import ParameterStructure
from stalk.pes.pes_function import NotEvaluatedException, PesFunction


class RelaxFunction(PesFunction):

    # Override the evaluate method to handle relaxation
    def _evaluate_structure(
        self,
        structure: ParameterSet,
        interactive: bool = False,
        reset_value=False,
        dep_jobs=[],
    ) -> None:
        """Evaluate the structure to be relaxed."""
        if reset_value:
            # Reset value
            structure.reset_value()
        else:
            # Try to load relaxed structure from disk
            params_relax = structure.load(path=structure.path, key='params_out')
            # Return if already evaluated or the cache loading succeeded
            if params_relax is not None:
                value, error = structure.value, structure.error
                structure.params = params_relax
                structure.value, structure.error = value, error
                print(f'{structure.path} is already relaxed.')
                # Load value, error if available
                structure.try_load_value()
                return
            # end if
        # end if
        if interactive:
            self._prompt([structure])
        # end if
        try:
            # Unlike elsewhere, the structure, value and error can be set here directly and
            # bypass using loaders. The data is transferred via the structure
            res = self.func(structure, **self.args)
            if isinstance(res, ParameterSet):
                value, error = res.value, res.error
                structure.params = res.params
                structure.value = value
                structure.error = error
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
                structure.save_value()
            # end if
            structure.save(path=structure.path, params_out=structure.params)
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
