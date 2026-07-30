#!/usr/bin/env python3
'''Class for saving/loading line-search optimizer data'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from stalk.io.txt_data import TxtData
from stalk.ls.error_surface import ErrorSurface
from stalk.ls.target_linesearch import TargetLineSearch


class OptimizerData():
    # Data file for the error surface correlated Gs
    Gs_file: TxtData = None
    # Data file for the error surface
    es_file: TxtData = None
    # Data file for the error surface X-mesh
    es_xmesh_file: TxtData = None
    # Data file for the error surface Y-mesh
    es_ymesh_file: TxtData = None
    # Data file for the optimal window size
    W_opt_file: TxtData = None
    # Data file for the optimal target errorbar
    sigma_opt_file: TxtData = None

    def __init__(
        self,
        label: str = 'ls',
    ):
        self.label = label
        self.Gs_file = TxtData(f'{label}_Gs.dat')
        self.es_file = TxtData(f'{label}_es.dat')
        self.es_xmesh_file = TxtData(f'{label}_es_X.dat')
        self.es_ymesh_file = TxtData(f'{label}_es_Y.dat')
        self.W_opt_file = TxtData(f'{label}_W_opt.dat')
        self.sigma_opt_file = TxtData(f'{label}_sigma_opt.dat')
    # end def

    def save(
        self,
        tls: TargetLineSearch,
        path: str,
        overwrite: bool = True
    ) -> None:
        self.Gs_file.save_result(path, tls.Gs, overwrite=overwrite)
        self.es_file.save_result(path, tls.error_surface.E_mat, overwrite=overwrite)
        self.es_xmesh_file.save_result(path, tls.error_surface.X_mat, overwrite=overwrite)
        self.es_ymesh_file.save_result(path, tls.error_surface.Y_mat, overwrite=overwrite)
        if tls.W_opt is not None and tls.sigma_opt is not None:
            self.W_opt_file.save_result(path, tls.W_opt, overwrite=overwrite)
            self.sigma_opt_file.save_result(path, tls.sigma_opt, overwrite=overwrite)
        # end if
    # end def

    def load(
        self,
        path: str,
    ) -> TargetLineSearch:
        result = TargetLineSearch(fit_kind='pf3')
        result.W_opt = self.W_opt_file.load_result(path, None)
        result.sigma_opt = self.sigma_opt_file.load_result(path, None)
        try:
            # Must have complete error surface data and Gs to load the error surface
            result.Gs = self.Gs_file.load_result(path, None)
            error_surface = ErrorSurface()
            error_surface._E_mat = self.es_file.load_result(path)
            error_surface._X_mat = self.es_xmesh_file.load_result(path)
            error_surface._Y_mat = self.es_ymesh_file.load_result(path)
            result._error_surface = error_surface
        except FileNotFoundError:
            pass
        # end try
        return result
    # end def

# end class
