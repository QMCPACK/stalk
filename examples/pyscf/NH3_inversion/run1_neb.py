#!/usr/bin/env python3

from ase.mep import NEB
from ase.optimize import BFGS

from stalk import interpolate_params
from stalk import XyzGeometry

from params import neb_image, pes_pbe
from run0_relax import structure_a, structure_b


# Generate base directory for NEB
basedir = 'neb'

# number of intermediate images
n_images = 3

traj_init = interpolate_params(structure_a, structure_b, n_images)
images = [neb_image(structure, pes_args=pes_pbe.args) for structure in traj_init]

# Try to load from disk
xyz = XyzGeometry(suffix='structure.xyz')
try:
    traj_neb = []
    for i in range(n_images + 2):
        res = xyz.load(f'{basedir}/image{i}/')
        image = structure_a.copy(pos=res.pos, label=f'image{i}')
        traj_neb.append(image)
    # end for
except FileNotFoundError:
    neb = NEB(images, climb=True)
    opt = BFGS(neb)
    opt.run(fmax=0.01)
    positions = opt.atoms.get_positions().reshape(-1, *structure_a.pos.shape)
    traj_neb = []
    for i, pyscf_image in enumerate(opt.atoms.images):
        image = structure_a.copy(pos=pyscf_image.positions, label=f'image{i}')
        traj_neb.append(image)
        xyz.write(image, f'{basedir}/image{i}/')
    # end for
# end try
