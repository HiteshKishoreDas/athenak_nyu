"""Periodic Favre filtering and SGS closure comparisons on regular synthetic grids.

This prepares a priori analysis. Synthetic checks do not establish closure accuracy
against a validated turbulent reference. No real reference dataset is bundled.
"""

import argparse
import json

import numpy as np


def box_filter(values, width):
    """Positive normalized top-hat filter over the first three periodic axes.

    Offsets start at -floor(width/2) and end at width-floor(width/2)-1. Even
    widths place the filtered field at a half-cell shifted center. All fields
    use the same offsets.
    Filtering preserves the fine grid sampling; it does not downsample.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.ndim < 3 or any(n == 0 for n in values.shape[:3]):
        raise ValueError("values must have three nonempty spatial axes")
    if (not isinstance(width, (int, np.integer)) or isinstance(width, bool)
            or width < 1 or width > min(values.shape[:3])):
        raise ValueError("width must be an integer in [1, smallest grid dimension]")
    if not np.isfinite(values).all():
        raise ValueError("filter input must be finite")
    result = values
    offsets = range(-(width // 2), width - width // 2)
    for axis in range(3):
        filtered = np.zeros_like(result)
        for offset in offsets:
            filtered += np.roll(result, -offset, axis=axis)
        result = filtered / width
    return result


def gaussian_kernel(n, width):
    """Positive periodic Gaussian; Delta_G = sqrt(12)*sigma in fine cells.

    This explicit second-moment convention is not certified paper-identical.
    Narrow sampled kernels differ from their continuous Fourier response.
    """
    if (not isinstance(n, (int, np.integer)) or n < 2
            or not np.isfinite(width) or not 0 < width <= n):
        raise ValueError("require n>=2 and 0<width<=n")
    sigma = width / np.sqrt(12)
    distance = (np.arange(n) + n // 2) % n - n // 2
    images = int(np.ceil(8 * sigma / n)) + 1
    weights = sum(np.exp(-0.5 * ((distance + image * n) / sigma)**2)
                  for image in range(-images, images + 1))
    return weights / weights.sum()


def gaussian_filter(values, width):
    """Same normalized positive separable convolution for every moment."""
    values = np.asarray(values, dtype=np.float64)
    if (values.ndim < 3 or len(set(values.shape[:3])) != 1
            or not np.isfinite(values).all()):
        raise ValueError("finite cubic spatial field required")
    response = np.fft.fft(gaussian_kernel(values.shape[0], width)).real
    transfer = (response[:, None, None] * response[None, :, None]
                * response[None, None, :])
    transfer = transfer.reshape(transfer.shape + (1,) * (values.ndim - 3))
    transformed = np.fft.fftn(values, axes=(0, 1, 2)) * transfer
    return np.fft.ifftn(transformed, axes=(0, 1, 2)).real


def velocity_gradient(velocity, dx):
    """A[...,i,j] = derivative_j velocity_i using periodic centered differences."""
    if not np.isfinite(dx) or dx <= 0:
        raise ValueError("dx must be finite and positive")
    a = np.empty(velocity.shape[:-1] + (3, 3), dtype=np.float64)
    for j in range(3):
        a[..., :, j] = (np.roll(velocity, -1, axis=j)
                       - np.roll(velocity, 1, axis=j)) / (2 * dx)
    return a


def favre_reference(rho, velocity, width, dx, filter_kind="box"):
    """Return filtered rho/u, exact tau/K and discrete exact production.

    Density is strictly positive; Favre division uses no denominator epsilon.
    Exact stress follows the negative Reynolds-stress convention of the paper.
    """
    rho = np.asarray(rho, dtype=np.float64)
    velocity = np.asarray(velocity, dtype=np.float64)
    if (rho.ndim != 3 or len(set(rho.shape)) != 1
            or velocity.shape != rho.shape + (3,)):
        raise ValueError("require rho shape (n,n,n) and velocity shape (n,n,n,3)")
    if (not np.isfinite(rho).all() or not np.isfinite(velocity).all()
            or np.any(rho <= 0)):
        raise ValueError("require finite velocity and finite strictly positive rho")
    filters = {"box": box_filter, "gaussian": gaussian_filter}
    if filter_kind not in filters:
        raise ValueError("filter_kind must be box or gaussian")
    spatial_filter = filters[filter_kind]
    rhobar = spatial_filter(rho, width)
    utilde = spatial_filter(rho[..., None] * velocity, width) / rhobar[..., None]
    tau = np.empty(rho.shape + (3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(i, 3):
            moment = spatial_filter(rho * velocity[..., i] * velocity[..., j], width)
            tau[..., i, j] = rhobar * utilde[..., i] * utilde[..., j] - moment
            tau[..., j, i] = tau[..., i, j]
    K = -0.5 * np.trace(tau, axis1=-2, axis2=-1)
    a = velocity_gradient(utilde, dx)
    # tau is symmetric, so tau:A = tau:S without storing a second tensor.
    production = np.einsum("...ij,...ij->...", tau, a)
    return {"rho": rhobar, "velocity": utilde, "tau": tau, "K": K,
            "gradient": a, "production": production}


def mixed_stress(rho, K, delta, a, c1=0.02, c2=0.7):
    """Independent Python implementation of Schmidt-Federrath Eq. (8).

    Exact zero gradient uses isotropic tau=-2K/3 I. Exact zero K gives zero.
    Input must be admissible; negative K is rejected rather than silently clipped.
    The straightforward tensor norm is intended for conditioned reference data.
    """
    if (not np.isfinite(delta) or delta <= 0
            or not np.isfinite(c1) or c1 < 0
            or not np.isfinite(c2) or not 0 <= c2 <= 1):
        raise ValueError("invalid Delta or mixed-closure coefficients")
    if (a.shape != rho.shape + (3, 3) or K.shape != rho.shape
            or not np.isfinite(a).all() or not np.isfinite(rho).all()
            or not np.isfinite(K).all() or np.any(rho <= 0) or np.any(K < 0)):
        raise ValueError("inadmissible mixed-closure state")
    q = np.einsum("...ij,...ij->...", a, a)
    if not np.isfinite(q).all() or np.any((q == 0) & np.any(a != 0, axis=(-2, -1))):
        raise ValueError("gradient norm range error; reference data outside safe range")
    d = np.trace(a, axis1=-2, axis2=-1)
    eddy = 2 * c1 * delta * np.sqrt(2) * np.sqrt(rho) * np.sqrt(K)
    tau = np.empty_like(a)
    nonzero = q > 0
    for i in range(3):
        for j in range(i, 3):
            strain = 0.5 * (a[..., i, j] + a[..., j, i])
            if i == j:
                strain = strain - d / 3
            aat = np.sum(a[..., i, :] * a[..., j, :], axis=-1)
            normalized = np.zeros_like(q)
            np.divide(aat, q, out=normalized, where=nonzero)
            value = eddy * strain - 2 * c2 * K * normalized
            if i == j:
                value -= (2 / 3) * (1 - c2) * K
                value = np.where(nonzero, value, -(2 / 3) * K)
            tau[..., i, j] = value
            tau[..., j, i] = value
    if not np.isfinite(tau).all():
        raise ValueError("mixed stress overflow; reference data outside safe range")
    return tau


def error_metrics(exact, model):
    """Signed model-minus-exact errors; undefined zero-reference ratios are null."""
    error = model - exact
    exact_rms = float(np.sqrt(np.mean(exact * exact)))
    rms_error = float(np.sqrt(np.mean(error * error)))
    return {"exact_mean": float(np.mean(exact)), "model_mean": float(np.mean(model)),
            "mean_signed_error": float(np.mean(error)),
            "rms_error": rms_error, "max_abs_error": float(np.max(np.abs(error))),
            "relative_rms_error": rms_error / exact_rms if exact_rms > 0 else None}


def production_metrics(exact, model, density):
    """Signed production statistics with explicit zero and density conventions.

    Correlation is Pearson's signed correlation, null for a constant field.
    Sign error excludes exact zero production but counts a zero prediction as
    incorrect for a nonzero reference. Density subsets split at the median;
    values equal to the median belong to the lower subset.
    """
    exact, model, density = (np.asarray(x, dtype=np.float64)
                             for x in (exact, model, density))
    if (exact.shape != model.shape or density.shape != exact.shape
            or exact.size == 0 or not np.isfinite(exact).all()
            or not np.isfinite(model).all() or not np.isfinite(density).all()
            or np.any(density <= 0)):
        raise ValueError("require same nonempty shapes and finite positive density")
    result = error_metrics(exact, model)
    x = exact - exact.mean()
    y = model - model.mean()
    denominator = np.sqrt(np.mean(x*x)) * np.sqrt(np.mean(y*y))
    nonzero = exact != 0
    result.update({"sample_count": int(exact.size),
                   "signed_correlation": float(np.mean(x*y)/denominator)
                   if denominator > 0 else None,
                   "exact_forward_fraction": float(np.mean(exact > 0)),
                   "exact_backscatter_fraction": float(np.mean(exact < 0)),
                   "model_forward_fraction": float(np.mean(model > 0)),
                   "model_backscatter_fraction": float(np.mean(model < 0)),
                   "false_sign_fraction_nonzero_exact": float(np.mean(
                       np.sign(model[nonzero]) != np.sign(exact[nonzero])))
                   if nonzero.any() else None})
    median = float(np.median(density))
    subsets = {}
    for name, mask in (("forward", exact > 0), ("backscatter", exact < 0),
                       ("density_at_or_below_median", density <= median),
                       ("density_above_median", density > median)):
        count = int(mask.sum())
        subsets[name] = {"sample_count": count,
                         "density_mean": float(density[mask].mean()) if count else None,
                         "errors": error_metrics(exact[mask], model[mask])
                         if count else None}
    result.update({"density_median": median, "subsets": subsets})
    return result


def compare_closure(rho, velocity, width, dx):
    reference = favre_reference(rho, velocity, width, dx)
    # Roundoff can make variance slightly negative. Do not hide this correction:
    # reject this comparison and require a better-conditioned reference calculation.
    if np.any(reference["K"] < 0):
        raise ValueError("negative filtered variance; no clipping applied")
    modeled_tau = mixed_stress(reference["rho"], reference["K"], width * dx,
                               reference["gradient"])
    modeled_production = np.einsum("...ij,...ij->...", modeled_tau,
                                    reference["gradient"])
    return {"stress": error_metrics(reference["tau"], modeled_tau),
            "production": production_metrics(reference["production"],
                                               modeled_production, reference["rho"]),
            "K_min": float(reference["K"].min()),
            "K_mean": float(reference["K"].mean()),
            "trace_identity_max_error": float(np.max(np.abs(
                np.trace(modeled_tau, axis1=-2, axis2=-1) + 2 * reference["K"])))}


def synthetic_metrics(n=16, width=3):
    """Small deterministic preparation sample, explicitly not a turbulent reference."""
    if not isinstance(n, int) or not 4 <= n <= 32:
        raise ValueError("synthetic grid must have 4 <= n <= 32")
    phase = 2 * np.pi * (np.arange(n) + 0.5) / n
    x, y, z = np.meshgrid(phase, phase, phase, indexing="ij")
    rho = 1 + 0.2 * np.cos(x) * np.cos(y)
    velocity = np.empty(rho.shape + (3,))
    velocity[..., 0] = np.sin(y) + 0.2 * np.cos(z)
    velocity[..., 1] = np.sin(z) + 0.2 * np.cos(x)
    velocity[..., 2] = np.sin(x) + 0.2 * np.cos(y)
    result = compare_closure(rho, velocity, width, 1 / n)
    result.update({"dataset": "deterministic synthetic periodic trigonometric fields",
                   "validated_turbulent_reference": False,
                   "claim": "filter and comparison pipeline preparation only",
                   "n": n, "filter_width_cells": width, "dx": 1 / n,
                   "Delta": width / n, "sampling": "filtered fine grid, no downsampling"})
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--width", type=int, default=3)
    args = parser.parse_args()
    print(json.dumps(synthetic_metrics(args.n, args.width), indent=2, allow_nan=False))
