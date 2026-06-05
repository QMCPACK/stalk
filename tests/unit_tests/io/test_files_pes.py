#!/usr/bin/env python

from numpy import savetxt

from stalk.io.files_pes import FilesPes
from stalk.io.pes_loader import PesLoader
from stalk.io.xyz_geometry import XyzGeometry
from stalk.io.files_pes import write_xyz_sigma
from stalk.util.util import match_to_tol
from tests.unit_tests.assets.h2o import get_structure_H2O

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_FilesPes(tmp_path):

    # Test default init
    pes = FilesPes()
    assert pes.func is write_xyz_sigma
    assert pes.args == {}
    assert pes.loader.suffix == 'energy.dat'
    assert isinstance(pes.loader, PesLoader)

    # Test evaluate
    s = get_structure_H2O()
    path = tmp_path / 'test0'
    # This creates the structure file but does not load the energy yet
    pes.evaluate(s, path=path)
    assert s.value is None
    assert s.error == 0.0
    res = XyzGeometry(suffix='structure.xyz').load(path)
    sigma_path = path / 'sigma.dat'
    assert not sigma_path.exists()
    assert match_to_tol(s.pos, res.get_pos())
    for e, e_ref in zip(s.elem, res.get_elem()):
        assert e == e_ref
    # end for
    assert s.value is None
    assert s.error == 0.0

    # Next, add energy
    value_ref = 1.0
    error_ref = 0.1
    savetxt(path / 'energy.dat', [value_ref, error_ref])
    pes.evaluate(s, path=path)
    assert match_to_tol(s.value, value_ref)
    assert match_to_tol(s.error, error_ref)

    # Test evaluate_all
    s2a = s.copy(label='2a')
    s2b = s.copy(label='2b')
    s2a.reset_value()
    s2b.reset_value()
    suffix = 'test_structure.xyz'
    sigma_suffix = 'test_sigma.dat'
    sigmas = [0.1, 0.2]
    pes.evaluate_all(
        [s2a, s2b],
        path=path,
        sigmas=sigmas,
        suffix=suffix,
        sigma_suffix=sigma_suffix,
    )
    values_ref = [1.1, 2.1]
    errors_ref = [0.11, 0.22]
    # Adding one structure energy but not the other
    savetxt(path / '2b/energy.dat', [values_ref[1], errors_ref[1]])
    pes.evaluate_all(
        [s2a, s2b],
        path=path,
        sigmas=sigmas,
        suffix=suffix,
        sigma_suffix=sigma_suffix
    )
    assert match_to_tol(s2b.value, values_ref[1])
    assert match_to_tol(s2b.error, errors_ref[1])
    assert s2a.value is None
    # Adding the remaining structure energy
    savetxt(path / '2a/energy.dat', [values_ref[0], errors_ref[0]])
    pes.evaluate_all(
        [s2a, s2b],
        path=path,
        sigmas=sigmas,
        suffix=suffix,
        sigma_suffix=sigma_suffix
    )
    assert match_to_tol(s2a.value, values_ref[0])
    assert match_to_tol(s2a.error, errors_ref[0])
    assert match_to_tol(s2b.value, values_ref[1])
    assert match_to_tol(s2b.error, errors_ref[1])

# end def
