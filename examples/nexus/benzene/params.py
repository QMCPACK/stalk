#!/usr/bin/env python3

from numpy import array, ndarray, sin, cos, pi, diag

from nexus import generate_pyscf, generate_qmcpack, job, obj
from nexus import generate_physical_system, generate_pw2qmcpack, generate_pwscf
from stalk.io.PesLoader import PesLoader
from stalk.io.util import load_energy
from structure import Structure

from stalk.params.util import mean_distances
from stalk.util.util import Bohr
from stalk.io.XyzGeometry import XyzGeometry
from stalk.nexus.NexusGeometry import NexusGeometry
from stalk.nexus.NexusPes import NexusPes
from stalk.nexus.QmcPes import QmcPes
from stalk.params.PesFunction import PesFunction
from stalk.util import EffectiveVariance

# This requires the following job arguments to be defined in local nxs.py
from nxs import pyscfjob, optjob, dmcjob, pwscfjob, p2qjob

# Pseudos (execute download_pseudos.sh in the working directory)
qmcpseudos = ['C.ccECP.xml', 'H.ccECP.xml']
scfpseudos = ['C.ccECP.upf', 'H.ccECP.upf']


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos: ndarray):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    # for easier comprehension, list particular atoms
    C0 = pos[0]
    C1 = pos[1]
    C2 = pos[2]
    C3 = pos[3]
    C4 = pos[4]
    C5 = pos[5]
    H0 = pos[6]
    H1 = pos[7]
    H2 = pos[8]
    H3 = pos[9]
    H4 = pos[10]
    H5 = pos[11]

    # for redundancy, calculate mean bond lengths
    # 0) from neighboring C-atoms
    r_CC = mean_distances([
        (C0, C1),
        (C1, C2),
        (C2, C3),
        (C3, C4),
        (C4, C5),
        (C5, C0)
    ])
    # 1) from corresponding H-atoms
    r_CH = mean_distances([
        (C0, H0),
        (C1, H1),
        (C2, H2),
        (C3, H3),
        (C4, H4),
        (C5, H5)
    ])
    params = array([r_CC, r_CH])
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params: ndarray):
    r_CC = params[0]
    r_CH = params[1]
    # place atoms on a hexagon in the xy-directions
    hex_xy = array([[cos(3 * pi / 6), sin(3 * pi / 6), 0.],
                    [cos(5 * pi / 6), sin(5 * pi / 6), 0.],
                    [cos(7 * pi / 6), sin(7 * pi / 6), 0.],
                    [cos(9 * pi / 6), sin(9 * pi / 6), 0.],
                    [cos(11 * pi / 6), sin(11 * pi / 6), 0.],
                    [cos(13 * pi / 6), sin(13 * pi / 6), 0.]])
    # C-atoms are one C-C length apart from origin
    pos_C = hex_xy * r_CC
    # H-atoms one C-H length apart from C-atoms
    pos_H = hex_xy * (r_CC + r_CH)
    pos = array([pos_C, pos_H])
    return pos
# end def


# return a 1-item list of Nexus jobs: SCF relaxation
def scf_relax_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        C=4,
        H=1
    )
    relax = generate_pyscf(
        template='../pyscf_relax.py',
        system=system,
        identifier='relax',
        job=job(**pyscfjob),
        path=path,
        mole=obj(
            verbose=4,
            ecp='ccecp',
            basis='ccpvdz',
            symmetry=False,
            cart=True
        ),
    )
    return [relax]
# end def


relax_pyscf = NexusGeometry(
    scf_relax_job,
    loader=XyzGeometry({'suffix': 'relax.xyz', 'c_pos': 1.0 / Bohr})
)


# Let us define an SCF PES job that is consistent with the earlier relaxation
def scf_pes_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        C=4,
        H=1,
    )
    scf = generate_pyscf(
        template='../pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path,
        mole=obj(
            verbose=4,
            ecp='ccecp',
            basis='ccpvdz',
            symmetry=False,
            cart=True
        ),
    )
    return [scf]
# end def


# Hessian based on the structural mappings
pes_pyscf = NexusPes(
    func=PesFunction(scf_pes_job),
    loader=PesLoader(load_energy, args={'suffix': 'energy.dat'})
)


# 4-5) Stochastic: Line-search
# return a simple 4-item DMC workflow for Nexus:
#   1: run SCF with norm-conserving ECPs to get orbitals
#   2: convert for QMCPACK
#   3: Optimize 2-body Jastrow coefficients
#   4: Run DMC with enough steps/block to meet the target errorbar sigma
#     it is important to first characterize the DMC run into var_eff
def dmc_pes_job(
    structure: Structure,
    path,
    sigma=None,
    samples=10,
    var_eff=None,
    dep_jobs=[],
    **kwargs
):
    # Estimate the relative number of samples needed
    if isinstance(var_eff, EffectiveVariance):
        dmcsteps = var_eff.get_samples(sigma)
    else:
        dmcsteps = samples
    # end if

    # For QMCPACK, use plane-waves for better performance
    axes = array([20., 20., 10.])
    structure.set_axes(diag(axes))
    structure.pos += axes / 2
    structure.kpoints = array([[0, 0, 0]])
    system = generate_physical_system(
        structure=structure,
        C=4,
        H=1,
    )
    scf = generate_pwscf(
        system=system,
        job=job(**pwscfjob),
        path=path + 'scf',
        pseudos=scfpseudos,
        identifier='scf',
        calculation='scf',
        input_type='generic',
        input_dft='pbe',
        nosym=False,
        nogamma=True,
        conv_thr=1e-9,
        mixing_beta=.7,
        ecutwfc=300,
        ecutrho=600,
        occupations='smearing',
        smearing='gaussian',
        degauss=0.0001,
        electron_maxstep=1000,
        kgrid=(1, 1, 1),  # needed to run plane-waves with Nexus
        kshift=(0, 0, 0,),
    )
    p2q = generate_pw2qmcpack(
        identifier='p2q',
        path=path + 'scf',
        job=job(**p2qjob),
        dependencies=[(scf, 'orbitals')],
    )
    system.bconds = 'nnn'
    reuse_jastrow = len(dep_jobs) > 0
    if reuse_jastrow:
        opt_cycles = 3
        minwalkers = 0.5
    else:
        opt_cycles = 8
        minwalkers = 0.1
    # end if
    opt = generate_qmcpack(
        system=system,
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
        J1_rcut=6.0,
        J2_size=8,
        J2_rcut=8.0,
        minmethod='oneshift',
        blocks=200,
        substeps=2,
        samples=100000,
        minwalkers=minwalkers,
    )
    if reuse_jastrow:
        opt.depends(dep_jobs[2])
    # end if
    dmc = generate_qmcpack(
        system=system,
        path=path + 'dmc',
        job=job(**dmcjob),
        dependencies=[(p2q, 'orbitals'), (opt, 'jastrow')],
        steps=dmcsteps,
        identifier='dmc',
        qmc='dmc',
        input_type='basic',
        pseudos=qmcpseudos,
        jastrows=[],
        walkers_per_rank=128,
        blocks=200,
        timestep=0.01,
        ntimesteps=1,
    )
    # Store the relative samples for printout
    dmc.samples = dmcsteps
    return [scf, p2q, opt, dmc]
# end def


# Configure a job generator and loader for the DMC PES
# -the suffix points to the correct nexus analyzer
# -the qmc_idx points to the correct QMC series (0: VMC; 1: DMC)
pes_dmc = NexusPes(
    dmc_pes_job,
    loader=QmcPes({'suffix': '/dmc/dmc.in.xml', 'qmc_idx': 1})
)
