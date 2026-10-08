"""Cell-centered closure diagnostic for aligned, block-filtered 3D snapshots.

Uses exact filtered K, so reservoir evolution cannot obscure the stress test.
This coarse-sampled tau:gradient diagnostic is not the solver's face production.
Native snapshot axes are (z,y,x); gradient columns always mean physical (x,y,z).
"""

import numpy as np

from filter_reference import error_metrics, mixed_stress, production_metrics
from filter_reference import velocity_gradient


def assess(data, coarse, K):
    """Return JSON-ready metrics and coarse-volume arrays for independent plots."""
    n = K.shape[0]
    if (K.shape != (n, n, n) or coarse.shape[:3] != K.shape
            or data.ndim != 4 or len(set(data.shape[:3])) != 1
            or data.shape[0] < n or data.shape[0] % n):
        raise ValueError("require aligned cubic data/coarse grids and matching K")
    r = data.shape[0] // n
    rho = data[..., 0]
    velocity = data[..., 1:4] / rho[..., None]
    favre = coarse[..., 1:4] / coarse[..., :1]
    expanded = favre.repeat(r, 0).repeat(r, 1).repeat(r, 2)
    fluctuation = velocity - expanded
    exact_tau = np.empty(K.shape + (3, 3))
    for i in range(3):
        for j in range(i, 3):
            moment = rho * fluctuation[..., i] * fluctuation[..., j]
            value = -moment.reshape(n, r, n, r, n, r).mean(axis=(1, 3, 5))
            exact_tau[..., i, j] = value
            exact_tau[..., j, i] = value
    # Existing reference helper expects spatial axes in physical x/y/z order.
    gradient = velocity_gradient(favre.transpose(2, 1, 0, 3), 1 / n)
    gradient = gradient.transpose(2, 1, 0, 3, 4)
    modeled_tau = mixed_stress(coarse[..., 0], K, 1 / n, gradient)
    exact_production = np.einsum("...ij,...ij->...", exact_tau, gradient)
    modeled_production = np.einsum("...ij,...ij->...", modeled_tau, gradient)
    metrics = dict(
        convention="negative Reynolds stress; positive tau:A feeds SGS K",
        sampling="aligned block filter; coarse centered gradient; whole volume",
        solver_face_production_reproduction=False,
        supplied_K="exact filtered kinetic variance",
        Delta=1 / n,
        filter_width_cells=r,
        stress=error_metrics(exact_tau, modeled_tau),
        deviatoric_stress=error_metrics(
            exact_tau + (2 / 3) * K[..., None, None] * np.eye(3),
            modeled_tau + (2 / 3) * K[..., None, None] * np.eye(3),
        ),
        production=production_metrics(
            exact_production, modeled_production, coarse[..., 0]
        ),
        exact_trace_K_max_error=float(np.max(abs(
            np.trace(exact_tau, axis1=-2, axis2=-1) + 2 * K))),
        model_trace_K_max_error=float(np.max(abs(
            np.trace(modeled_tau, axis1=-2, axis2=-1) + 2 * K))),
    )
    return dict(metrics=metrics, fields=dict(
        exact_tau=exact_tau, modeled_tau=modeled_tau, gradient=gradient,
        exact_production=exact_production, modeled_production=modeled_production,
    ))


def self_test():
    data = np.zeros((2, 2, 2, 7))
    data[..., 0] = [1, 3]
    data[..., 1] = [1, 9]
    coarse = data.mean(axis=(0, 1, 2), keepdims=True)
    result = assess(data, coarse, np.full((1, 1, 1), 0.75))
    assert abs(result["fields"]["exact_tau"][0, 0, 0, 0, 0] + 1.5) < 1e-14
    assert result["metrics"]["exact_trace_K_max_error"] < 1e-14
    assert result["metrics"]["model_trace_K_max_error"] < 1e-14
    assert result["metrics"]["production"]["rms_error"] == 0
    n = 4
    data = np.zeros((8, 8, 8, 7))
    data[..., 0] = 1
    phase = 2 * np.pi * (np.arange(8) + 0.5) / 8
    data[..., 1] = np.sin(phase)[None, :, None]
    coarse = data.reshape(n, 2, n, 2, n, 2, 7).mean(axis=(1, 3, 5))
    expanded = coarse[..., 1].repeat(2, 0).repeat(2, 1).repeat(2, 2)
    K = (0.5 * (data[..., 1] - expanded)**2).reshape(
        n, 2, n, 2, n, 2).mean(axis=(1, 3, 5))
    result = assess(data, coarse, K)
    gradient = result["fields"]["gradient"]
    expected = n * np.sin(2 * np.pi / n) * np.cos(np.pi / 8) * np.cos(
        2 * np.pi * (np.arange(n) + 0.5) / n)
    assert np.max(abs(gradient[..., 0, 1] - expected[None, :, None])) < 1e-13
    gradient[..., 0, 1] = 0
    assert np.max(abs(gradient)) == 0
    assert result["metrics"]["exact_trace_K_max_error"] < 1e-14
    print("Aligned covariance, K trace and physical gradient-axis checks pass")


if __name__ == "__main__":
    self_test()
