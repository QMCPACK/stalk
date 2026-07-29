#!/usr/bin/env python3

from numpy import ndarray

from nexus import generate_qmcpack, job, obj
from nexus import generate_physical_system, generate_pw2qmcpack, generate_pwscf

from stalk.util import Bohr
from stalk import Parameter
from stalk import PwscfGeometry
from stalk import PwscfPes
from stalk import NexusGeometry
from stalk import NexusPes
from stalk import QmcPes
from stalk import NexusStructure

# This requires the following job arguments to be defined in local nxs.py
# Copy examples/nexus/nxs_template.py to ./nxs.py and edit accordingly
from nxs import pwscfjob, optjob, dmcjob, p2qjob  # type: ignore

# Pseudos (execute download_pseudos.sh in the working directory)
softpseudos = ['C.pbe_v1.2.uspp.F.upf']
scfpseudos = ['C.ccECP.upf']
qmcpseudos = ['C.ccECP.xml']


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos: ndarray, axes: ndarray):
    # Redundancy helps in finding silly mistakes in the parameter mappings
    a = Parameter([
        axes[0, 0],
        axes[0, 1],
        axes[1, 1],
        axes[1, 2],
        axes[2, 0],
        axes[2, 2],
    ], label='a')
    return [a]
# end def


def backward(params: ndarray):
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
def scf_vcrelax_job(structure: NexusStructure, **kwargs):
    system = generate_physical_system(
        structure=structure.get_nexus_structure(),
        C=4,
    )
    relax = generate_pwscf(
        system=system,
        job=job(**pwscfjob),
        path=structure.path,
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
def scf_pes_job(structure: NexusStructure, **kwargs):
    system = generate_physical_system(
        structure=structure.get_nexus_structure(),
        C=4,
    )
    scf = generate_pwscf(
        system=system,
        job=job(**pwscfjob),
        path=structure.path,
        identifier='scf',
        calculation='scf',
        **scf_args
    )
    return [scf]
# end def


# Nexus generator for DMC PES workflow
def dmc_pes_job(
    structure: NexusStructure,
    dep_jobs=[],
    tile_opt=1,
    rcut=3.0,
    twist_grid=(1, 1, 1),
    kshift=(0, 0, 0),
    **kwargs
):
    path = structure.path
    dmcsteps = structure.samples if structure.samples is not None else 10
    structure_nx = structure.get_nexus_structure()
    # The original structure represents the primitive cell
    primcell = generate_physical_system(
        structure=structure_nx,
        C=4
    )
    # Create spherical tiling to mitigate finite-size effects
    tiled = structure_nx.tile_opt(tile_opt)
    tiled.add_symmetrized_kmesh(kgrid=(1, 1, 1), kshift=kshift)
    supercell = generate_physical_system(structure=tiled, C=4)
    # Twisted supercell
    tiled_twisted = structure_nx.tile_opt(tile_opt)
    tiled_twisted.add_symmetrized_kmesh(kgrid=twist_grid, kshift=kshift)
    twisted_supercell = generate_physical_system(structure=tiled_twisted, C=4)

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
    nscft = generate_pwscf(
        system=twisted_supercell,
        job=job(**pwscfjob),
        path=path + 'nscft',
        identifier='nscf',
        calculation='nscf',
        wf_collect=True,
        nosym=True,
        nogamma=True,
        dependencies=[(scf, 'charge_density')],
        **scf_qmc_args
    )
    p2qt = generate_pw2qmcpack(
        job=job(**p2qjob),
        path=path + 'nscft',
        identifier='p2q',
        write_psir=False,
        dependencies=[(nscft, 'orbitals')],
    )
    dmc = generate_qmcpack(
        system=twisted_supercell,
        path=path + 'dmc',
        job=job(**dmcjob),
        dependencies=[(p2qt, 'orbitals'), (opt, 'jastrow')],
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
    )
    return [scf, nscf, p2q, opt, nscft, p2qt, dmc]
# end def


# Finally, wrap the Nexus job generators as defined above with appropriate loader arguments
# to be used as finalized relaxation/PES recipes.
relax_pwscf = NexusGeometry(
    scf_vcrelax_job,
    # Nexus PwscfAnalyzer returns PWSCF relaxed structure in Angstrom
    loader=PwscfGeometry(suffix='vcrelax.in', scale=Bohr**-1)
)
pes_pwscf = NexusPes(
    scf_pes_job,
    # Nexus PwscfAnalyzer returns PWSCF energy
    loader=PwscfPes(suffix='scf.in')
)
pes_dmc = NexusPes(
    dmc_pes_job,
    # Nexus QmcpackAnalyzer returns DMC energy for the first time-step after walker
    # generation, so at index->1
    loader=QmcPes(suffix='dmc/dmc.in', qmc_idx=1)
)
