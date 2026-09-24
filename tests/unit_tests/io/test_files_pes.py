#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import savetxt, loadtxt

from stalk.io.files_pes import FilesPes
from stalk.pes.pes_loader import PesLoader
from stalk.io.xyz_geometry import XyzGeometry
from stalk.io.files_pes import write_xyz_sigma
from stalk.util.util import match_to_tol
from tests.unit_tests.assets.h2o import get_structure_H2O


def test_FilesPes(tmp_path):

    # Test default init
    pes = FilesPes()
    assert pes.func is write_xyz_sigma
    assert pes.args == {}
    assert pes.loader.suffix == 'value.out'
    assert isinstance(pes.loader, PesLoader)

    # Test evaluate
    s = get_structure_H2O()
    sigma = 0.1
    s.sigma = sigma
    path = tmp_path / 'test0'
    # This creates the structure, sigma files but does not load the energy yet
    pes.evaluate(s, path=path)
    # Sigma should be created
    sigma_path = s.path / 'sigma.in'
    assert sigma_path.exists()
    assert loadtxt(sigma_path, ndmin=1)[0] == sigma
    # Geometry should be created
    structure_path = s.path / 'structure.xyz'
    res = XyzGeometry(suffix='structure.xyz').load(structure_path)
    assert match_to_tol(s.pos, res.get_pos())
    for e, e_ref in zip(s.elem, res.get_elem()):
        assert e == e_ref
    # end for
    # Params should also be created
    params_path = s.path / 'params.in'
    assert params_path.exists()
    # Value shoud not be available yet
    assert s.value is None
    assert s.error == 0.0
    # Next, write energies to disk and evaluate again to load the values
    value_path = s.path / 'value.out'
    error_path = s.path / 'error.out'
    value_ref = 1.0
    error_ref = 0.1
    savetxt(value_path, [value_ref])
    savetxt(error_path, [error_ref])
    pes.evaluate(s, path=path)
    assert match_to_tol(s.value, value_ref)
    assert match_to_tol(s.error, error_ref)

    # Test evaluate_all
    s2a = s.copy(label='2a')
    s2b = s.copy(label='2b')
    s2a.reset_value()
    s2b.reset_value()
    pes.evaluate_all(
        [s2a, s2b],
        path=path,
    )
    values_ref = [1.1, 2.1]
    errors_ref = [0.11, 0.22]
    # Adding one structure energy but not the other
    savetxt(path / '2b/value.out', [values_ref[1]])
    savetxt(path / '2b/error.out', [errors_ref[1]])
    pes.evaluate_all(
        [s2a, s2b],
        path=path,
    )
    assert match_to_tol(s2b.value, values_ref[1])
    assert match_to_tol(s2b.error, errors_ref[1])
    assert s2a.value is None
    # Adding the remaining structure energy
    savetxt(path / '2a/value.out', [values_ref[0]])
    savetxt(path / '2a/error.out', [errors_ref[0]])
    pes.evaluate_all(
        [s2a, s2b],
        path=path,
    )
    assert match_to_tol(s2a.value, values_ref[0])
    assert match_to_tol(s2a.error, errors_ref[0])
    assert match_to_tol(s2b.value, values_ref[1])
    assert match_to_tol(s2b.error, errors_ref[1])

# end def
