# Routine to rebin piecewise powerlaw distribution

import numpy as np
from scipy.optimize import brentq

from dust_consts import *
import utils as ut
import ic


def rebin(
        dist, 
        slopes, 
        centers, 
        edges, 
        shift, 
        flat_in_bin=False,
    ):

    under, over = 0.0, 0.0  # grains/cm^3

    # * Add ghost bins for over and underflow
    dist_pad = np.zeros(np.shape(dist)[0] + 2, dtype=float)  # 1/um
    slopes_pad = np.zeros(np.shape(dist)[0] + 2, dtype=float)  # 1/um
    bin_pad = np.zeros(np.shape(centers)[0] + 2, dtype=float)  # 1/um
    edges_pad = np.ones(np.shape(edges)[0] + 2, dtype=float)

    shifted_bin_pad = np.ones(np.shape(centers)[0] + 2, dtype=float)
    shifted_edges_pad = np.zeros(np.shape(edges)[0] + 2, dtype=float)

    dist_pad[1:-1] = dist
    slopes_pad[1:-1] = slopes
    bin_pad[1:-1] = centers  # um
    edges_pad[1:-1] = edges  # um

    bin_pad[0] = ic.a * 1e-2
    bin_pad[-1] = ic.b * 1e2

    edges_pad[0] = ic.a * 1e-4
    edges_pad[-1] = ic.b * 1e4

    shifted_bin_pad = bin_pad + shift  # um
    shifted_edges_pad = edges_pad + shift  # um

    bin_sizes = edges_pad[1:] - edges_pad[:-1]  # um

    def exact_integral(x1, x2, mass=False):

        integral_func = ut.bin_num_integral
        if mass:
            integral_func = ut.bin_mass_integral

        return integral_func(
            x1,
            x2,
            dist_pad,
            slopes_pad,
            shifted_edges_pad,
            shifted_bin_pad,
            amax=shifted_edges_pad[0],
            amin=shifted_edges_pad[-1],
        )

    new_dist = np.zeros_like(dist_pad)  # 1/um
    new_slopes = np.zeros_like(slopes_pad)  # 1/um
    # Solve only physical bins. i=0 and i=-1 are ghost bins.
    for i in range(1, len(dist_pad) - 1):
        left_edge = edges_pad[i]
        right_edge = edges_pad[i + 1]

        def integral_in_bin(p):
            return bin_pad[i] * ut.integrate_power(
                p,
                left_edge / bin_pad[i],
                right_edge / bin_pad[i],
            )

        num_in_bin = exact_integral(left_edge, right_edge)
        mass_in_bin = exact_integral(left_edge, right_edge, mass=True)

        if not flat_in_bin:

            if num_in_bin==0.0:
                p = 0
            else:
                def solve_lhs(p):
                    ac3_bin = mass_in_bin / num_in_bin / K_dust
                    ac3_bin /= bin_pad[i]**3
                    return (integral_in_bin(p + 3) / integral_in_bin(p)) - ac3_bin

                print(f"{solve_lhs(-30) = }")
                print(f"{solve_lhs(30) = }")

                p = brentq(solve_lhs, -30, 30)

            new_slopes[i] = p

            # print(f"Solved p: {p}")

        new_dist[i] = num_in_bin / integral_in_bin(new_slopes[i])

    under = exact_integral(shifted_edges_pad[0], edges_pad[1], mass=True)
    over = exact_integral(edges_pad[-2], shifted_edges_pad[-1], mass=True)

    return new_dist[1:-1], new_slopes[1:-1], under, over


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    Nbins = 10

    test_dist_ind = -2

    # Demo test: apply a small positive shift and visualize the rebinned result.
    edges = np.logspace(-3, -1, num=Nbins + 1)
    bin = np.sqrt(edges[1:] * edges[:-1])
    dist = np.ones_like(bin) * 1e12 * (bin / bin[0]) ** test_dist_ind
    shift = 3.0e-3

    slopes = np.zeros_like(dist) + test_dist_ind

    # rebin() uses module-level a and b to define ghost-bin bounds.
    a_test = edges[0]
    b_test = edges[-1]

    rebinned, slopes_rebin, under, over = rebin(
        dist=dist, slopes=slopes, centers=bin, edges=edges, shift=shift,
        flat_in_bin=False,
    )

    # Total particles before and after rebinning
    def bin_sum(dist, slopes, shift=0, mass=False):
        if mass:
            int_func = ut.bin_mass_integral
        else:
            int_func = ut.bin_num_integral

        return int_func(
            a_test+shift,
            b_test+shift,
            dist,
            slopes,
            edges+shift,
            bin+shift,
            amax=a_test+shift,
            amin=b_test+shift,
        )


    tot_num_before = bin_sum(dist, slopes, shift=shift)
    tot_num_after = bin_sum(rebinned, slopes_rebin)
    tot_mass_before = bin_sum(dist, slopes, mass=True, shift=shift)
    tot_mass_after = bin_sum(rebinned, slopes_rebin, mass=True)


    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    ut.plot_piecewise_powerlaw_distribution(
        dust_bin_edges=edges,
        dust_bin=bin,
        dust_dist=dist,
        slope_bin=slopes,
        info_str="Original distribution",
        ax=axs[0],
        show=False,
    )
    ut.plot_piecewise_powerlaw_distribution(
        dust_bin_edges=edges,
        dust_bin=bin,
        dust_dist=rebinned,
        slope_bin=slopes_rebin,
        info_str="Rebinned distribution",
        ax=axs[1],
        show=False,
    )
    fig.suptitle(f"rebin.py demo test (shift = {shift:.2e} um)")
    fig.tight_layout()
    plt.show()
    plt.close()

    print("rebin.py demo rebinning complete.")
    print(f"underflow mass: {under:.6e} g/cm^3")
    print(f"overflow mass:  {over:.6e} g/cm^3\n")

    print("=============================")
    print(f"{tot_num_before = }")
    print(f"{tot_num_after = }\n")
    print(f"{tot_mass_before = }")
    print(f"{tot_mass_after = }\n")

    print("=============================")
    print(f"Mass conservation check:")
    print(f"{tot_mass_before = }")
    print(f"{under+over+tot_mass_after = }\n")

    residual = under+over+tot_mass_after -  tot_mass_before
    print(f"{residual = }")
