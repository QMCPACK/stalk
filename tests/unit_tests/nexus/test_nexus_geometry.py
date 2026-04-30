#!/usr/bin/env python3

from pytest import raises, warns

from stalk.nexus.nexus_geometry import NexusGeometry
from stalk.nexus.nexus_structure import NexusStructure
from stalk.util.util import match_to_tol

from ..assets.test_jobs import nxs_generic_pes, TestGeometryLoader
from ..assets.h2o import pos_H2O, elem_H2O
from ..assets.diamond import pos_diamond, elem_diamond, axes_diamond


__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_NexusGeometry(tmp_path):

    # Test empty (should fail)
    with raises(TypeError):
        NexusGeometry()
    # end with

    pos_H2O_shifted = pos_H2O * 1.2
    s = NexusStructure(
        label='label',
        pos=pos_H2O_shifted,
        elem=elem_H2O,
        units='A'
    )
    # Test successful relaxation job, non-periodic
    pes = NexusGeometry(
        nxs_generic_pes,
        args={'pes_variable': 'relax_h2o'},
        loader=TestGeometryLoader(),
    )
    assert not s.analyzed
    assert not s.generated
    pes.relax(s, path=str(tmp_path) + "/relax")
    # The structure should have been updated to eqm values
    assert match_to_tol(s.pos, pos_H2O)
    # Test unsuccessful relaxation job
    s.set_position(pos_H2O_shifted)
    pes.loader.args['produce_fail'] = True
    with warns(UserWarning):
        pes.relax(s, path=str(tmp_path) + "/relax_fail")
        # The position should be unchanged
        assert match_to_tol(s.pos, pos_H2O_shifted)
    # end with

    # Test successful relaxation, periodic
    pos_diamond_shifted = pos_diamond * 1.1
    axes_diamond_shifted = axes_diamond * 1.1
    s2 = NexusStructure(
        label='label',
        pos=pos_diamond_shifted,
        axes=axes_diamond_shifted,
        elem=elem_diamond,
        units='A'
    )
    pes2 = NexusGeometry(
        nxs_generic_pes,
        args={'pes_variable': 'relax_diamond'},
        loader=TestGeometryLoader(produce_fail=False)
    )
    pes2.relax(s2, path=str(tmp_path) + "/relax_periodic")
    # The sctructure should have been updated to eqm values
    assert match_to_tol(s2.pos, pos_diamond)
    assert match_to_tol(s2.axes, axes_diamond)
    # Test unsuccessful relaxation, periodic
    pes2.loader.args['produce_fail'] = True
    s2.set_position(pos_diamond_shifted, axes_diamond_shifted)
    with warns(UserWarning):
        pes2.relax(s2, path=str(tmp_path) + "/relax_periodic_fail")
        # The positions and axes should be unchanged
        assert match_to_tol(s2.pos, pos_diamond_shifted)
        assert match_to_tol(s2.axes, axes_diamond_shifted)
    # end with

# end def
