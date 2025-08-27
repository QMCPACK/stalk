#!/usr/bin/env python3

from numpy import ndarray

from nexus import generate_qmcpack, job, obj
from nexus import generate_physical_system, generate_pw2qmcpack, generate_pwscf
from structure import Structure

from stalk.nexus.PwscfGeometry import PwscfGeometry
from stalk.nexus.PwscfPes import PwscfPes
from stalk.util.util import Bohr
from stalk.nexus.NexusGeometry import NexusGeometry
from stalk.nexus.NexusPes import NexusPes
from stalk.nexus.QmcPes import QmcPes
from stalk.util import EffectiveVariance

# This requires the following job arguments to be defined in local nxs.py
from nxs import pwscfjob, optjob, dmcjob, p2qjob

# Pseudos (execute download_pseudos.sh in the working directory)
softpseudos = ['C.pbe_v1.2.uspp.F.upf']
scfpseudos = ['C.ccECP.upf']
qmcpseudos = ['C.ccECP.xml']


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos: ndarray, axes: ndarray):
    from stalk.params.util import mean_param
    from numpy import array
    # Redundancy helps in finding silly mistakes in the parameter mappings
    a = mean_param([
        axes[0, 0],
        axes[0, 1],
        axes[1, 1],
        axes[1, 2],
        axes[2, 0],
        axes[2, 2],
    ])
    return array([a])
# end def


def backward(params):
    a = params[0]
    axes = [
        [a, a, 0],
        [0, a, a],
        [a, 0, a]
    ]
    pos = [
        [0.0, 0.0, 0.0],
        [a / 2, a / 2, a / 2]
    ]
    return pos, axes
# end def


# Define common PWSCF arguments to keep consistent between relaxation and PES.
scf_args = obj(
    pseudos=softpseudos,
    ecutwfc=80,
    ecutrho=300,
    kgrid=(8, 8, 8),
    kshift=(0, 0, 0,),
    input_type='generic',
    input_dft='pbe',
    occupations='smearing',
    smearing='gaussian',
    conv_thr=1e-9,
    degauss=0.0001,
    nosym=False,
    mixing_beta=.7,
    electron_maxstep=200,
)


# Nexus generator for PWSCF relaxation workflow
def scf_vcrelax_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        C=4,
    )
    relax = generate_pwscf(
        system=system,
        job=job(**pwscfjob),
        path=path,
        identifier='vcrelax',
        calculation='vc-relax',
        forc_conv_thr=1e-4,
        ion_dynamics='bfgs',
        press=0.0,
        press_conv_thr=0.4,
        **scf_args
    )
    return [relax]
# end def


# Nexus generator for SCF PES workflow
def scf_pes_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        C=4,
    )
    scf = generate_pwscf(
        system=system,
        job=job(**pwscfjob),
        path=path,
        identifier='scf',
        calculation='scf',
        **scf_args
    )
    return [scf]
# end def


# Nexus generator for DMC PES workflow
def dmc_pes_job(
    structure: Structure,
    path,
    sigma=None,
    samples=10,
    var_eff=None,
    tile_opt=1,
    dep_jobs=[],
    rcut=3.0,
    **kwargs
):
    # Estimate the relative number of samples needed
    if isinstance(var_eff, EffectiveVariance):
        dmcsteps = var_eff.get_samples(sigma)
    else:
        dmcsteps = samples
    # end if

    # The original structure represents the primitive cell
    primcell = generate_physical_system(
        structure=structure,
        C=4
    )
    # Create spherical tiling to mitigate finite-size effects
    tiled = structure.tile_opt(tile_opt)
    # In this example, we only consider twistnumber=0
    tiled.add_symmetrized_kmesh(kgrid=(1, 1, 1), kshift=(0, 0, 0))
    supercell = generate_physical_system(
        structure=tiled,
        C=4,
    )

    # Check if reusing Jastrows
    reuse_jastrow = len(dep_jobs) > 0 and tiled.rwigner(1) > rcut
    if reuse_jastrow:
        opt_cycles = 4
        minwalkers = 0.5
    else:
        opt_cycles = 10
        minwalkers = 0.2
        rcut = None
    # end if

    # Common SCF args
    scf_qmc_args = obj(
        pseudos=scfpseudos,
        input_type='generic',
        input_dft='pbe',
        conv_thr=1e-9,
        mixing_beta=.7,
        ecutwfc=300,
        ecutrho=600,
        occupations='smearing',
        smearing='gaussian',
        degauss=0.0001,
        electron_maxstep=200,
    )
    scf = generate_pwscf(
        system=primcell,
        job=job(**pwscfjob),
        path=path + 'scf',
        identifier='scf',
        calculation='scf',
        kgrid=(8, 8, 8),
        kshift=(0, 0, 0,),
        wf_collect=False,
        **scf_qmc_args
    )
    nscf = generate_pwscf(
        system=supercell,
        job=job(**pwscfjob),
        path=path + 'nscf',
        identifier='nscf',
        calculation='nscf',
        wf_collect=True,
        nosym=True,
        nogamma=True,
        dependencies=[(scf, 'charge_density')],
        **scf_qmc_args
    )
    p2q = generate_pw2qmcpack(
        identifier='p2q',
        path=path + 'nscf',
        job=job(**p2qjob),
        write_psir=False,
        dependencies=[(nscf, 'orbitals')],
    )
    opt = generate_qmcpack(
        system=supercell,
        path=path + 'opt',
        job=job(**optjob),
        dependencies=[(p2q, 'orbitals')],
        cycles=opt_cycles,
        identifier='opt',
        qmc='opt',
        input_type='basic',
        pseudos=qmcpseudos,
        J2=True,
        J1_size=6,
        J1_rcut=rcut,
        J2_size=8,
        J2_rcut=rcut,
        minmethod='oneshift',
        blocks=200,
        substeps=2,
        samples=20000,
        meshfactor=0.8,
        minwalkers=minwalkers,
        # Using the first twist
        twistnum=0,
    )
    if reuse_jastrow:
        opt.depends(dep_jobs[3], 'jastrow')
    # end if
    dmc = generate_qmcpack(
        system=supercell,
        path=path + 'dmc',
        job=job(**dmcjob),
        dependencies=[(p2q, 'orbitals'), (opt, 'jastrow')],
        steps=dmcsteps,
        identifier='dmc',
        qmc='dmc',
        input_type='basic',
        pseudos=qmcpseudos,
        jastrows=[],
        meshfactor=0.8,
        walkers_per_rank=32,
        blocks=200,
        timestep=0.01,
        ntimesteps=1,
        nonlocalmoves=True,
        twistnum=0,
    )
    # Store the relative samples for printout
    dmc.samples = dmcsteps
    return [scf, nscf, p2q, opt, dmc]
# end def


# Finally, wrap the Nexus job generators as defined above with appropriate loader arguments
# to be used as finalized relaxation/PES recipes.
relax_pwscf = NexusGeometry(
    scf_vcrelax_job,
    # Nexus PwscfAnalyzer returns PWSCF relaxed structure in Angstrom
    loader=PwscfGeometry({'suffix': 'vcrelax.in', 'c_pos': Bohr})
)
pes_pwscf = NexusPes(
    scf_pes_job,
    # Nexus PwscfAnalyzer returns PWSCF energy
    loader=PwscfPes({'suffix': 'scf.in'})
)
pes_dmc = NexusPes(
    dmc_pes_job,
    # Nexus QmcpackAnalyzer returns DMC energy for the first time-step after walker
    # generation, so at index->1
    loader=QmcPes({'suffix': '/dmc/dmc.in.xml', 'qmc_idx': 1})
)
