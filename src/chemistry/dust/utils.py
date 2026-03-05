import numpy as np
from dust_consts import *


def integrate_power(slope, x0, x1, a_c=1, num_c=1, tol=1e-8):
    # Return the integral of (x**slope) between x0 and x1
    eps = slope + 1.0

    # near a = -1 -> log
    if abs(eps) < tol:
        return np.log(x1 / x0)

    # stable difference using expm1
    # factor out exp(eps*log(x0)) to avoid cancellation when x1 ~ x0
    t = eps * np.log(x1 / x0)
    return np.exp(eps * np.log(x0)) * np.expm1(t) / eps

def integrate_linear(slope, x0, x1, a_c, num_c):
    # Return the integral of (1+ slope*a_c/center*(x-1))

    # For number integrals, n(a)=n_c + slope*(a-a_c), so the dimensionless
    # coefficient in x=a/a_c space is slope*a_c/n_c.
    result =  0.0 if num_c==0 else slope*a_c/num_c
    result *= 0.5*(x1**2 - x0**2) - x1 + x0
    result += x1-x0
    
    return result

def integrate_linear_mass(slope, x0, x1, a_c, num_c):
    # Return the integral of (1+ slope*a_c/center*(x-1))
    result =  0.0 if num_c==0 else slope*(a_c**4)/num_c
    result *= 0.2*(x1**5 - x0**5) - 0.25*(x1**4-x0**4) 
    result += 0.25*(x1**4-x0**4)

    return result


def bin_num_integral(a1, a2, dist, slope, edges, centers, amax, amin, type="linear"):

    integrate_dict = {}
    integrate_dict["linear"] = integrate_linear
    integrate_dict["powerlaw"] = integrate_power
    integrate_dict["linear_mass"] = integrate_linear_mass

    integrate_fn = integrate_dict[type]


    ai = max(a1, amin)
    af = min(a2, amax)
    if af <= ai:
        return 0.0

    ind_i = np.searchsorted(edges, ai) - 1
    ind_f = np.searchsorted(edges, af) - 1
    # ind_i = np.clip(ind_i, 0, len(dist) - 1)
    # ind_f = np.clip(ind_f, 0, len(dist) - 1)

    integral = 0.0

    if ind_i == ind_f:
        integral += dist[ind_i] * centers[ind_i] * integrate_fn(
            slope[ind_i],
            ai / centers[ind_i],
            af / centers[ind_f],
            a_c = centers[ind_i],
            num_c = dist[ind_i],
        )
        return integral

    # Add the left edge
    integral += dist[ind_i] * centers[ind_i] * integrate_fn(
        slope[ind_i],
        ai / centers[ind_i],
        edges[ind_i + 1] / centers[ind_i],
        a_c = centers[ind_i],
        num_c = dist[ind_i],
    )

    # Add the right edge
    integral += dist[ind_f] * centers[ind_f] * integrate_fn(
        slope[ind_f],
        edges[ind_f] / centers[ind_f],
        af / centers[ind_f],
        a_c = centers[ind_f],
        num_c = dist[ind_f],
    )

    # integrate the bins in between
    if ind_f > ind_i + 1:
        integral += np.sum(
            dist[ind_i + 1 : ind_f]
            * centers[ind_i + 1 : ind_f]
            * np.vectorize(integrate_fn)(
                slope[ind_i + 1 : ind_f],
                edges[ind_i + 1 : ind_f] / centers[ind_i + 1 : ind_f],
                edges[ind_i + 2 : ind_f + 1] / centers[ind_i + 1 : ind_f],
                a_c = centers[ind_i + 1 : ind_f],
                num_c = dist[ind_i + 1 : ind_f],
            )
        )

    return integral


def bin_mass_integral(a1, a2, dist, slope, edges, centers, amax, amin, type="linear"):

    new_slope = slope.copy()

    if type=="powerlaw":
        new_slope += 3
    elif type=="linear":
        type = "linear_mass"

    return (
        K_dust 
        * bin_num_integral(
            a1,
            a2,
            dist * centers**3,
            new_slope,
            edges,
            centers,
            amax,
            amin,
            type=type,
        )
    )

def plot_piecewise_powerlaw_distribution(
    dust_bin_edges,
    dust_bin,
    dust_dist,
    slope_bin,
    info_str="",
    ax=None,
    samples_per_bin=20,
    show=True,
):
    import matplotlib.pyplot as plt

    n_bins = len(dust_bin)
    if len(dust_bin_edges) != n_bins + 1:
        raise ValueError("dust_bin_edges must have length len(dust_bin) + 1")
    if len(dust_dist) != n_bins or len(slope_bin) != n_bins:
        raise ValueError("dust_dist and slope_bin must match len(dust_bin)")

    created_ax = ax is None
    if created_ax:
        plt.figure()
        ax = plt.gca()

    for i in range(n_bins):
        a1, a2 = dust_bin_edges[i], dust_bin_edges[i + 1]
        a_center = dust_bin[i]
        n_center = dust_dist[i]
        slope = slope_bin[i]

        a_segment = np.logspace(np.log10(a1), np.log10(a2), samples_per_bin)
        n_segment = n_center * (a_segment / a_center) ** slope

        ax.plot(a_segment, n_segment, linestyle="dashed", color=f"C{i % 10}")
        ax.plot(
            a_segment,
            n_center * np.ones_like(a_segment),
            linestyle="solid",
            color="C0",
        )
        ax.axvline(a1, color="C1", linestyle="dotted")

    ax.axvline(dust_bin_edges[-1], color="C1", linestyle="dotted", label="Bin edges")
    ax.set_xlabel("Grain size (micrometers)")
    ax.set_ylabel("Number density")
    ax.set_title(info_str)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid()
    ax.legend()

    if created_ax:
        ax.figure.suptitle("Initial dust size distribution per cm$^3$")
        if show:
            plt.show()
    return ax


def plot_piecewise_linear_distribution(
    dust_bin_edges,
    dust_bin,
    dust_dist,
    slope_bin,
    info_str="",
    ax=None,
    samples_per_bin=20,
    show=True,
):
    import matplotlib.pyplot as plt

    n_bins = len(dust_bin)
    if len(dust_bin_edges) != n_bins + 1:
        raise ValueError("dust_bin_edges must have length len(dust_bin) + 1")
    if len(dust_dist) != n_bins or len(slope_bin) != n_bins:
        raise ValueError("dust_dist and slope_bin must match len(dust_bin)")

    created_ax = ax is None
    if created_ax:
        plt.figure()
        ax = plt.gca()

    for i in range(n_bins):
        a1, a2 = dust_bin_edges[i], dust_bin_edges[i + 1]
        a_center = dust_bin[i]
        n_center = dust_dist[i]
        slope = slope_bin[i]

        a_segment = np.logspace(np.log10(a1), np.log10(a2), samples_per_bin)
        n_segment = n_center + slope * (a_segment - a_center)

        ax.plot(a_segment, n_segment, linestyle="dashed", color=f"C{i % 10}")
        ax.plot(
            a_segment,
            n_center * np.ones_like(a_segment),
            linestyle="solid",
            color="C0",
        )
        ax.axvline(a1, color="C1", linestyle="dotted")

    ax.axvline(dust_bin_edges[-1], color="C1", linestyle="dotted", label="Bin edges")
    ax.set_xlabel("Grain size (micrometers)")
    ax.set_ylabel("Number density")
    ax.set_title(info_str)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid()
    ax.legend()

    if created_ax:
        ax.figure.suptitle("Initial dust size distribution per cm$^3$")
        if show:
            plt.show()
    return ax


if __name__ == "__main__":
    # Test 1: near expo = -1 should take the logarithmic branch.
    x0 = 0.2
    x1 = 1.3
    got = integrate_power(-1.0 + 1e-12, x0, x1)
    expected = np.log(x1 / x0)
    print("[test1] integrate_power near expo=-1")
    print(f"[test1] x0={x0}, x1={x1}")
    print(f"[test1] expected={expected:.16e}, got={got:.16e}")
    assert np.isclose(got, expected, rtol=1e-12, atol=0.0)
    print("[test1] PASS")

    # Test 2: single-bin integration should match direct analytic value.
    edges = np.array([1.0, 2.0])
    centers = np.array([1.5])
    dist = np.array([3.0])
    slope = np.array([2.0])
    a1 = 1.2
    a2 = 1.8
    amax = 0.0
    amin = 10.0

    left = a1 / centers[0]
    right = a2 / centers[0]
    expected_num = dist[0] * centers[0] * (right**3 - left**3) / 3.0
    got_num = bin_num_integral(a1, a2, dist, slope, edges, centers, amax, amin)
    print("[test2] bin_num_integral single-bin analytic check")
    print(
        f"[test2] a1={a1}, a2={a2}, dist={dist[0]}, slope={slope[0]}, "
        f"center={centers[0]}"
    )
    print(f"[test2] expected_num={expected_num:.16e}, got_num={got_num:.16e}")
    assert np.isclose(got_num, expected_num, rtol=1e-12, atol=0.0)
    print("[test2] PASS")

    got_mass = bin_mass_integral(a1, a2, dist, slope, edges, centers, amax, amin)
    q = slope[0]
    expected_mass = (
        K_dust
        * dist[0]
        * centers[0] ** 4
        * (right ** (q + 4.0) - left ** (q + 4.0))
        / (q + 4.0)
    )
    print("[test3] bin_mass_integral single-bin analytic check")
    print(f"[test3] expected_mass={expected_mass:.16e}, got_mass={got_mass:.16e}")
    assert np.isclose(got_mass, expected_mass, rtol=1e-12, atol=0.0)
    print("[test3] PASS")

    # Test 4: multi-bin integration should match piecewise analytic sum.
    edges_mb = np.array([1.0, 2.0, 4.0, 8.0])
    centers_mb = np.array([1.5, 3.0, 6.0])
    dist_mb = np.array([2.0, 1.0, 0.5])
    slope_mb = np.array([0.0, -0.5, 1.0])
    a1_mb = 1.4
    a2_mb = 6.3
    amax_mb = 0.0
    amin_mb = 10.0

    def primitive(expo, x):
        eps = expo + 1.0
        if abs(eps) < 1e-14:
            return np.log(x)
        return x**eps / eps

    expected_mb = 0.0
    for i in range(len(dist_mb)):
        seg_l = max(a1_mb, edges_mb[i])
        seg_r = min(a2_mb, edges_mb[i + 1])
        if seg_r <= seg_l:
            continue
        u_l = seg_l / centers_mb[i]
        u_r = seg_r / centers_mb[i]
        expected_mb += (
            dist_mb[i]
            * centers_mb[i]
            * (primitive(slope_mb[i], u_r) - primitive(slope_mb[i], u_l))
        )

    got_mb = bin_num_integral(
        a1_mb, a2_mb, dist_mb, slope_mb, edges_mb, centers_mb, amax_mb, amin_mb
    )
    print("[test4] bin_num_integral multi-bin piecewise check")
    print(
        f"[test4] range=({a1_mb}, {a2_mb}), bins={len(dist_mb)}, "
        f"expected_mb={expected_mb:.16e}, got_mb={got_mb:.16e}"
    )
    assert np.isclose(got_mb, expected_mb, rtol=1e-12, atol=0.0)
    print("[test4] PASS")

    got_mb_mass = bin_mass_integral(
        a1_mb, a2_mb, dist_mb, slope_mb, edges_mb, centers_mb, amax_mb, amin_mb
    )

    expected_mb_mass = 0.0
    for i in range(len(dist_mb)):
        seg_l = max(a1_mb, edges_mb[i])
        seg_r = min(a2_mb, edges_mb[i + 1])
        if seg_r <= seg_l:
            continue
        u_l = seg_l / centers_mb[i]
        u_r = seg_r / centers_mb[i]
        expected_mb_mass += (
            K_dust
            * dist_mb[i]
            * centers_mb[i] ** 4
            * (primitive(slope_mb[i] + 3.0, u_r) - primitive(slope_mb[i] + 3.0, u_l))
        )
    print("[test5] bin_mass_integral multi-bin analytic check")
    print(
        f"[test5] expected_mb_mass={expected_mb_mass:.16e}, "
        f"got_mb_mass={got_mb_mass:.16e}"
    )
    assert np.isclose(got_mb_mass, expected_mb_mass, rtol=1e-12, atol=0.0)
    print("[test5] PASS")

    print("utils.py self-tests passed")
