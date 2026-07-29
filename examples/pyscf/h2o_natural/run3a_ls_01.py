#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches

from stalk import LineSearchIteration

from params import pes_dict, co_dict
from run2_surrogate import surrogates


# Run line-searches between PES combinations
lsis = {}
for xc_srg, pes_srg in pes_dict.items():
    lsis[xc_srg] = {}
    for xc_ls, pes_ls in pes_dict.items():
        path = f'ls_01/{xc_srg}-{xc_ls}'
        structure = surrogates[xc_srg].structure.copy()
        if xc_srg == xc_ls:
            structure.shift_params([0.1, -0.2])
        # end if
        lsi = LineSearchIteration(
            surrogate=surrogates[xc_srg],
            structure=structure,
            path=path,
        )
        for i in range(4):
            lsi.propagate(pes_ls, i, add_sigma=True)
        # end for
        # Evaluate the latest eqm structure
        lsi[-1].evaluate_eqm(pes_ls, add_sigma=True)
        print(f'Line-search ({xc_ls} + noise) on {xc_srg} surrogate:')
        print(lsi)
        print(surrogates[xc_ls].structure.params)
        print('^^Reference params^^')
        lsis[xc_srg][xc_ls] = lsi
    # end for
# end for


# Plot
if __name__ == '__main__':
    to_deg = 180 / np.pi
    for xc_ls, pes_ls in pes_dict.items():
        f, ax = plt.subplots()
        # Plot reference parameters
        co = co_dict[xc_ls]
        params = surrogates[xc_ls].structure.params
        epsilon_p = surrogates[xc_ls].epsilon_p
        ellipse = patches.Ellipse(
            (params[0], params[1] * to_deg),
            2 * epsilon_p[0],
            2 * epsilon_p[1] * to_deg,
            color=co,
            alpha=0.2
        )
        ax.add_patch(ellipse)
        ax.plot(params[0], params[1] * to_deg, marker='o', color=co, linestyle='None')
        # Plot ls trajectories
        for xc_srg, pes_srg in pes_dict.items():
            lsi = lsis[xc_srg][xc_ls]
            params = [lsi.pls(0).structure.params]
            params_err = [lsi.pls(0).structure.params_err]
            for pls in lsi.pls_list:
                if pls.evaluated:
                    params.append(pls.structure_next.params)
                    params_err.append(pls.structure_next.params_err)
                # end if
            # end for
            params = np.array(params)
            params_err = np.array(params_err)
            ax.errorbar(
                params[:, 0],
                params[:, 1] * to_deg,
                xerr=params_err[:, 0],
                yerr=params_err[:, 1] * to_deg,
                marker='.',
                linestyle=':',
                color=co_dict[xc_srg],
                label=f'{xc_ls} on {xc_srg}'
            )
        # end for
        plt.legend()
    # end for
    plt.xlabel('Bond length (A)')
    plt.ylabel('Bond angle (deg)')
    plt.show()
# end if
