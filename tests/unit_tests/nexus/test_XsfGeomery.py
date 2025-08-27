#!/usr/bin/env python3

from numpy import array
from pytest import raises

from stalk.params import GeometryResult
from stalk.nexus.XsfGeometry import XsfGeometry
from stalk.nexus.NexusStructure import NexusStructure
from stalk.params.ParameterStructure import ParameterStructure
from stalk.util.util import match_to_tol
from ..assets.gese import pos_GeSe, axes_GeSe, elem_GeSe

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"


def test_XsfGeometry(tmp_path):

    # default args: (suffix: structure.xsf)
    loader = XsfGeometry()

    # Test failing to load a file (default suffix: structure.xsf)
    with raises(FileNotFoundError):
        loader.load('tests/unit_tests/assets/pwscf_relax')
    # end with

    # Test loading a reference file
    res = loader.load('tests/unit_tests/assets/pwscf_relax', suffix='relax.xsf')

    # These are hardcoded in 'tests/unit_tests/assets/pwscf_relax/relax.xsf'
    elem_ref = ['V', 'Se', 'Se']
    pos_ref = array([
        [-0.00000000, 1.92228769, 9.98849984],
        [1.66474997, 0.96114385, 11.56470960],
        [1.66474997, 0.96114385, 8.41229008]
    ])
    axes_ref = array([
        [3.32949995, 0.00000000, 0.00000000],
        [-1.66474997, 2.88343154, 0.00000000],
        [0.00000000, 0.00000000, 19.97699968]
    ])

    assert isinstance(res, GeometryResult)
    assert match_to_tol(res.pos, pos_ref)
    assert match_to_tol(res.axes, axes_ref)
    assert all(elem_ref == res.elem)

    # Test writer
    writer = XsfGeometry()

    with raises(TypeError):
        writer.write(ParameterStructure())
    # end with
    # XSF will be written in Angstrom units
    s = NexusStructure(pos=pos_GeSe, axes=axes_GeSe, elem=elem_GeSe, units='A')
    writer.write(s, tmp_path, suffix='testfile.xsf')

    # Check by loading
    write_res = writer.load(tmp_path, suffix='testfile.xsf')
    assert isinstance(write_res, GeometryResult)
    assert match_to_tol(write_res.pos, pos_GeSe)
    assert match_to_tol(write_res.axes, axes_GeSe)
    assert all(elem_GeSe == write_res.elem)

# end def
