import numpy as np  

from dust_consts import *

#!!!! Make it conserve mass !!!!
def rebin(dist, bin, edges, shift):

    under, over = 0.0, 0.0  # grains/cm^3
    # shift_back = -shift


    # * Add ghost bins for over and underflow
    dist_pad = np.zeros(np.shape(dist)[0] + 2, dtype=float)  # 1/um
    shifted_bin_pad = np.ones(np.shape(bin)[0] + 2, dtype=float)
    edges_pad = np.ones(np.shape(bin)[0] + 2, dtype=float)
    shifted_edges_pad = np.zeros(np.shape(edges)[0] + 2, dtype=float)

    dist_pad[1:-1] = dist
    shifted_edges_pad[1:-1] = edges + shift  # um
    shifted_bin_pad[1:-1] = bin + shift  # um

    shifted_edges_pad[0] = a * 1e-4 + shift
    shifted_edges_pad[-1] = 1e6 * b + shift
    shifted_bin_pad[0] = a * 1e-2 + shift
    shifted_bin_pad[-1] = b * 1e3 + shift

    edges_pad = shifted_edges_pad - shift  # um
    bin_sizes = edges_pad[1:] - edges_pad[:-1]  # um

    def exact_integral(x1, x2, mass=False):
        # Integrate piecewise-constant dist_pad on shifted bins.
        if x2 <= x1:
            return 0.0

        i1 = np.searchsorted(shifted_edges_pad, x1, side="right") - 1
        i2 = np.searchsorted(shifted_edges_pad, x2, side="left") - 1

        i1 = np.clip(i1, 0, len(dist_pad) - 1)
        i2 = np.clip(i2, 0, len(dist_pad) - 1)

        if i1 == i2:
            return dist_pad[i1] * (x2 - x1)

        next_left_edge = shifted_edges_pad[i1 + 1]
        prev_right_edge = shifted_edges_pad[i2]

        mult=np.ones_like(shifted_bin_pad)  # dimensionless
        if mass:
            mult = 4 * np.pi * rho_gr_um * (shifted_bin_pad**3) / 3  # g/um^3

        integral = 0.0
        integral += dist_pad[i1] * (next_left_edge - x1) * mult[i1]
        integral += dist_pad[i2] * (x2 - prev_right_edge) * mult[i2]
        if i2 > i1 + 1:
            integral += np.sum(dist_pad[i1 + 1 : i2] * bin_sizes[i1 + 1 : i2] * mult[i1 + 1 : i2])
        return integral

    shift_dist = np.zeros_like(dist_pad)  # 1/um
    for i in range(1, len(edges_pad) - 1):
        left_edge = edges_pad[i]
        right_edge = edges_pad[i + 1]
        shift_dist[i] = exact_integral(left_edge, right_edge) / bin_sizes[i]

    under = exact_integral(shifted_edges_pad[0], edges_pad[1], mass=True)
    over = exact_integral(edges_pad[-2], shifted_edges_pad[-1], mass=True)

    return shift_dist[1:-1], under, over
