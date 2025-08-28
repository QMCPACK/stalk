#!/usr/bin/env python3

from numpy import array

from nexus import generate_pyscf, generate_qmcpack, job, obj
from nexus import generate_physical_system, generate_convert4qmc
from stalk.io.util import load_energy
from structure import Structure

from stalk.params.util import distance
from stalk.io.XyzGeometry import XyzGeometry
from stalk.nexus.NexusGeometry import NexusGeometry
from stalk.nexus.NexusPes import NexusPes
from stalk.nexus.QmcPes import QmcPes
from stalk.util import EffectiveVariance

# This requires the following job arguments to be defined in local nxs.py
from nxs import pyscfjob, optjob, dmcjob

# Pseudos (execute download_pseudos.sh in the working directory)
qmcpseudos = ['C.ccECP.xml']


# Forward mapping: produce parameter values from an array of atomic positions
def forward(pos):
    pos = pos.reshape(-1, 3)  # make sure of the shape
    C0 = pos[0]
    C1 = pos[1]
    d = distance(C0, C1)
    params = array([d])
    return params
# end def


# Backward mapping: produce array of atomic positions from parameters
def backward(params, **kwargs):
    d = params[0]
    C0 = [0.0, 0.0, -d / 2]
    C1 = [0.0, 0.0, +d / 2]
    pos = array([C0, C1]).flatten()
    return pos
# end def


# Define common SCF Mole arguments to keep consistent between relaxation and PES.
scf_mole_args = obj(
    spin=4,  # For SCF convergence, important to study eqm spin state
    verbose=4,
    ecp='ccecp',
    basis='ccpvdz',
    symmetry=False,
)


# Nexus generator for SCF relaxation workflow
def scf_relax_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        C=4
    )
    relax = generate_pyscf(
        template='../pyscf_relax.py',
        system=system,
        identifier='relax',
        job=job(**pyscfjob),
        path=path,
        mole=scf_mole_args
    )
    return [relax]
# end def


# Nexus generator for SCF PES workflow
def scf_pes_job(structure: Structure, path, **kwargs):
    system = generate_physical_system(
        structure=structure,
        C=4,
    )
    scf = generate_pyscf(
        template='../pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path,
        mole=scf_mole_args
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
    dep_jobs=[],
    **kwargs
):
    # Estimate the relative number of samples needed
    if isinstance(var_eff, EffectiveVariance):
        dmcsteps = var_eff.get_samples(sigma)
    else:
        dmcsteps = samples
    # end if

    # Check if reusing Jastrows
    reuse_jastrow = len(dep_jobs) > 0
    if reuse_jastrow:
        opt_cycles = 4
        minwalkers = 0.5
    else:
        opt_cycles = 10
        minwalkers = 0.2
    # end if

    # Center the structure for QMCPACK
    system = generate_physical_system(
        structure=structure,
        C=4,
        net_spin=4,
    )
    scf = generate_pyscf(
        template='../pyscf_pes.py',
        system=system,
        identifier='scf',
        job=job(**pyscfjob),
        path=path + 'scf',
        mole=obj(
            spin=4,
            verbose=4,
            ecp='ccecp',
            basis='ccecp-ccpvqz',  # Use larger basis to promote QMC performance
            symmetry=False,
        ),
        save_qmc=True,
    )
    c4q = generate_convert4qmc(
        identifier='c4q',
        path=path + 'scf',
        job=job(cores=1),
        dependencies=(scf, 'orbitals'),
    )
    opt = generate_qmcpack(
        system=system,
        path=path + 'opt',
        job=job(**optjob),
        dependencies=[(c4q, 'orbitals')],
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
        opt.depends(dep_jobs[2], 'jastrow')
    # end if
    dmc = generate_qmcpack(
        system=system,
        path=path + 'dmc',
        job=job(**dmcjob),
        dependencies=[(c4q, 'orbitals'), (opt, 'jastrow')],
        steps=dmcsteps,
        identifier='dmc',
        qmc='dmc',
        input_type='basic',
        pseudos=qmcpseudos,
        jastrows=[],
        walkers_per_rank=256,
        blocks=100,
        timestep=0.01,
        ntimesteps=1,
    )
    # Store the relative samples for printout
    dmc.samples = dmcsteps
    return [scf, c4q, opt, dmc]
# end def


# Finally, wrap the Nexus job generators as defined above with appropriate loader arguments
# to be used as finalized relaxation/PES recipes.
relax_pyscf = NexusGeometry(
    scf_relax_job,
    # pyscf_relax.py is configured to output relaxed geometry in relax.xyz
    loader=XyzGeometry({'suffix': 'relax.xyz'})
)
pes_pyscf = NexusPes(
    scf_pes_job,
    # pyscf_pes.py is configured to output SCF energy in energy.dat
    load_func=load_energy
)
pes_dmc = NexusPes(
    dmc_pes_job,
    # Nexus QmcpackAnalyzer returns DMC energy for the first time-step after walker
    # generation, so at index->1
    loader=QmcPes({'suffix': '/dmc/dmc.in.xml', 'qmc_idx': 1})
)
