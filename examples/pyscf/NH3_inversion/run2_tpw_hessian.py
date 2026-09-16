#!/usr/bin/env python3

from stalk import TransitionPathway

from params import pes_pbe
from run1_neb import traj_neb

tpw = TransitionPathway(
    path='tpw',
    images=traj_neb
)
tpw.calculate_hessians(pes=pes_pbe)

if __name__ == '__main__':
    for image in tpw.images:
        print(image.structure.params)
        print(image.hessian)
    # end for
# end if
