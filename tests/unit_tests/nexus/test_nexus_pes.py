#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pytest import raises

from stalk.params.pes_function import NotEvaluatedException
from stalk.util.util import match_to_tol
from stalk.nexus.nexus_pes import NexusPes
from stalk.nexus.nexus_structure import NexusStructure
from stalk import EffectiveVarianceMap

from ..assets.test_jobs import nxs_generic_pes, TestLoader
from ..assets.h2o import pes_H2O, pos_H2O, elem_H2O, forward_H2O


def test_NexusPes(tmp_path):

    s = NexusStructure(
        label='label',
        pos=pos_H2O,
        elem=elem_H2O,
        units='A'
    )

    # Test empty (should fail)
    with raises(TypeError):
        NexusPes()
    # end with

    # 1a: Test successful generation of jobs
    pes = NexusPes(
        nxs_generic_pes,
        args={'pes_variable': 'h2o'},
        loader=TestLoader()
    )
    pes.loader.args = {}  # ensure that other tests are not interfering with the args
    assert not pes.disable_failed
    assert not pes.bundle_jobs
    pes.evaluate(s, path=str(tmp_path) + '/nosigma', sigma=0.0)
    assert s.generated
    assert len(s.jobs) == 1
    assert s.finished
    E_original = pes_H2O(pos_H2O)[0]
    assert match_to_tol(s.value, E_original)
    assert match_to_tol(s.error, 0.0)
    # 1b: Test successful generation of jobs with noise
    sigma = 0.1
    s.reset_value()
    pes.evaluate(s, path=str(tmp_path) + '/sigma', sigma=sigma, add_sigma=True)
    # The value must have shifted
    assert not match_to_tol(s.value, E_original)
    assert match_to_tol(s.error, 0.1)
    # 1c: Test unsuccessful loading of jobs
    s.reset_value()
    pes.loader.args = {'produce_fail': True}  # make the loader fail
    with raises(NotEvaluatedException):
        pes.evaluate(s, path=str(tmp_path) + '/fail', sigma=0.0)
    # end with
    assert s.enabled
    assert match_to_tol(s.error, 0.0)
    # 1d: Test disabling of failed jobs
    s.reset_value()
    pes.disable_failed = True
    pes.evaluate(s, path=str(tmp_path) + '/disable_fail')
    assert not s.enabled

    # 2: Test evaluate all
    pes.loader.args = {}  # rest loader args
    sigmas = [0.1, 0.2, 0.3, 0.4]
    s2a = s.copy(pos=pos_H2O * 0.9, label='s2a')
    s2eqm1 = s.copy(pos=pos_H2O, label='eqm')
    s2b = s.copy(pos=pos_H2O * 1.1, label='s2b')
    s2eqm2 = s.copy(pos=pos_H2O, label='eqm')
    structures = [s2a, s2eqm1, s2b, s2eqm2]
    pes.evaluate_all(
        structures,
        path=str(tmp_path) + '/eval_all',
        sigmas=sigmas,
        add_sigma=True
    )
    # All are generated but the later duplicate eqm should match the first one
    assert all([s.generated for s in structures])
    assert structures[-1].jobs == []
    assert match_to_tol([s.error for s in structures], sigmas)
    # Re-evaluation makes no difference
    pes.evaluate_all(
        structures,
        path=str(tmp_path) + '/eval_all',
        sigmas=sigmas,
        add_sigma=True
    )
    # 2a: Test dep_jobs
    s2de = s.copy(pos=pos_H2O, label='eqm')
    s2dd = s.copy(pos=pos_H2O * 1.2, label='dep')
    s2dd.path = str(tmp_path) + '/dep'
    dep_jobs = nxs_generic_pes(s2dd)
    assert not dep_jobs[0].finished
    pes.evaluate_all(
        [s2de],
        dep_jobs=dep_jobs,
        path=str(tmp_path) + '/eval_dep'
    )
    assert s2de.generated
    assert not s2dd.generated
    assert dep_jobs[0].finished

    # 3: Test bundle
    sigmas = [0.1, 0.2]
    pes_bundle = NexusPes(
        nxs_generic_pes,
        args={'pes_variable': 'fail'},
        loader=TestLoader(),
        bundle_jobs=True
    )
    assert pes_bundle.bundle_jobs
    s1 = s.copy(pos=pos_H2O * 0.9, label='bundle1')
    s2 = s.copy(pos=pos_H2O * 1.1, label='bundle2')
    structures_bundle = [s1, s2]
    pes.evaluate_all(
        structures_bundle,
        path=str(tmp_path) + '/bundle',
        sigmas=sigmas,
        add_sigma=True
    )
    assert all([s.generated for s in structures_bundle])
    assert match_to_tol([s.error for s in structures_bundle], sigmas)

    # 4: Test effective variance map
    samples = 15
    pes = NexusPes(
        nxs_generic_pes,
        args={'pes_variable': 'evm'},
        loader=TestLoader()
    )
    pes.loader.args = {}  # ensure that other tests are not interfering with the args
    s_evm = s.copy()
    evm = pes.get_var_eff_map(
        structure=s_evm,
        path=str(tmp_path) + '/evm_test',
        samples=samples
    )
    assert isinstance(evm, EffectiveVarianceMap)
    assert len(evm) == 1
    assert evm.scaling_map[0][0] is s_evm
    assert match_to_tol(evm.scaling_map[0][1].var_eff, s_evm.error**2 * samples)

    s_evm0 = s_evm.copy(label='evm_test0', pos=pos_H2O * 0.9)
    sigma = 0.0001
    s_evm.forward = forward_H2O
    s_evm0.forward = forward_H2O
    samples_ref = evm.get_samples(s_evm0, error=sigma)
    pes.evaluate(
        s_evm0,
        path=str(tmp_path) + '/evm_test0',
        var_eff_map=evm,
        sigma=sigma,
    )
    assert s_evm0.samples == samples_ref
    assert len(evm) == 2
    assert evm.scaling_map[1][0] is s_evm0
    # TODO: The new estimated value is not externally validated
    error_new = 0.004
    samples_new = 277
    assert evm.get_samples(s_evm0, error=error_new) == samples_new
# end def
