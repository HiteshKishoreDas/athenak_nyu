"""Independent hand/Fourier/variance checks, frozen absolute tolerance 1e-12.

Coverage: F01-F05/F17 pipeline preparation. These synthetic checks do not prove
closure accuracy on a validated turbulent reference.
"""

import unittest

import numpy as np

from filter_reference import (box_filter, compare_closure, favre_reference,
                              error_metrics, mixed_stress, synthetic_metrics,
                              velocity_gradient, production_metrics)

TOL = 1e-12


class FilterReferenceTests(unittest.TestCase):
    def test_constant_preservation_and_positivity(self):
        values = np.full((4, 4, 4), 3.25)
        for width in range(1, 5):
            np.testing.assert_allclose(box_filter(values, width), values,
                                       atol=TOL, rtol=0)
        values[0, 0, 0] = 0
        self.assertGreaterEqual(box_filter(values, 3).min(), 0)
        self.assertLessEqual(box_filter(values, 3).max(), 3.25)

    def test_density_weighted_hand_case(self):
        # Equal-volume cells: rho=(1,3), ux=(1,3). rhobar=2, uFavre=2.5,
        # filter(rho ux^2)=14, so tau_xx=2*2.5^2-14=-1.5 and K=0.75.
        rho = np.ones((2, 2, 2))
        rho[1] = 3
        velocity = np.zeros((2, 2, 2, 3))
        velocity[0, ..., 0] = 1
        velocity[1, ..., 0] = 3
        ref = favre_reference(rho, velocity, 2, 0.5)
        np.testing.assert_allclose(ref["rho"], 2, atol=TOL, rtol=0)
        np.testing.assert_allclose(ref["velocity"][..., 0], 2.5, atol=TOL, rtol=0)
        np.testing.assert_allclose(ref["tau"][..., 0, 0], -1.5, atol=TOL, rtol=0)
        np.testing.assert_allclose(ref["K"], 0.75, atol=TOL, rtol=0)
        np.testing.assert_allclose(ref["production"], 0, atol=TOL, rtol=0)

    def test_sinusoid_transfer_function(self):
        n = 16
        phase = 2 * np.pi * np.arange(n) / n
        values = np.broadcast_to(np.sin(phase)[:, None, None], (n, n, n))
        gain = (1 + 2 * np.cos(2 * np.pi / n)) / 3
        np.testing.assert_allclose(box_filter(values, 3), gain * values, atol=TOL, rtol=0)

    def test_even_filter_center_and_derivative_commutation(self):
        n = 8
        phase = 2*np.pi*np.arange(n)/n
        values = np.broadcast_to(np.sin(phase)[:, None, None], (n, n, n))
        expected = np.broadcast_to((np.cos(np.pi/n)*np.sin(phase-np.pi/n))
                                    [:, None, None], (n, n, n))
        np.testing.assert_allclose(box_filter(values, 2), expected, atol=TOL, rtol=0)
        velocity = np.zeros((n, n, n, 3))
        velocity[..., 0] = values
        filtered_gradient = velocity_gradient(box_filter(velocity, 2), 1/n)
        gradient_filtered = box_filter(velocity_gradient(velocity, 1/n), 2)
        np.testing.assert_allclose(filtered_gradient, gradient_filtered,
                                   atol=TOL, rtol=0)

    def test_variance_trace_sign_and_galilean_invariance(self):
        n = 8
        phase = 2 * np.pi * (np.arange(n) + 0.5) / n
        x, y, z = np.meshgrid(phase, phase, phase, indexing="ij")
        rho = 1 + 0.2 * np.cos(x)
        velocity = np.stack((np.sin(x), np.cos(y), np.sin(z)), axis=-1)
        ref = favre_reference(rho, velocity, 3, 1 / n)
        boosted = favre_reference(rho, velocity + [5, -3, 2], 3, 1 / n)
        np.testing.assert_allclose(ref["tau"], boosted["tau"], atol=TOL, rtol=0)
        np.testing.assert_allclose(ref["K"], boosted["K"], atol=TOL, rtol=0)
        np.testing.assert_allclose(np.trace(ref["tau"], axis1=-2, axis2=-1),
                                   -2 * ref["K"], atol=TOL, rtol=0)
        self.assertGreaterEqual(ref["K"].min(), -TOL)
        self.assertLessEqual(np.linalg.eigvalsh(ref["tau"]).max(), TOL)

    def test_derivative_component_order_and_signed_error_metrics(self):
        n = 8
        phase = 2 * np.pi * np.arange(n) / n
        velocity = np.zeros((n, n, n, 3))
        velocity[..., 1] = np.sin(phase)[:, None, None]
        a = velocity_gradient(velocity, 1 / n)
        expected = np.broadcast_to((n * np.sin(2 * np.pi / n) * np.cos(phase))
                                    [:, None, None], (n, n, n))
        np.testing.assert_allclose(a[..., 1, 0], expected, atol=TOL, rtol=0)
        a[..., 1, 0] = 0
        np.testing.assert_allclose(a, 0, atol=TOL, rtol=0)
        metrics = error_metrics(np.array([-1., 2.]), np.array([0., 0.]))
        self.assertAlmostEqual(metrics["mean_signed_error"], -0.5, delta=TOL)
        self.assertAlmostEqual(metrics["rms_error"], np.sqrt(2.5), delta=TOL)

    def test_mixed_stress_backscatter_hand_case(self):
        # A=diag(.2,-.1,-.1), Q=.06 and S*=A. The nonlinear contribution
        # to tau:S is -.14; eddy contribution is 12*.02*sqrt(2)*.1^2.
        a = np.diag([0.2, -0.1, -0.1]).reshape((1, 1, 1, 3, 3))
        rho = np.ones((1, 1, 1))
        tau = mixed_stress(rho, rho, 1, a)
        expected_diag = (2 * 0.02 * np.sqrt(2) * np.array([0.2, -0.1, -0.1])
                         - 2 * 0.7 * np.array([4, 1, 1]) / 6 - 0.2)
        np.testing.assert_allclose(np.diag(tau[0, 0, 0]), expected_diag, atol=TOL, rtol=0)
        production = float(np.sum(tau * a))
        expected = 12 * 0.02 * np.sqrt(2) * 0.1**2 - 0.14
        self.assertAlmostEqual(production, expected, delta=TOL)
        self.assertLess(production, 0)

    def test_zero_gradient_zero_energy_and_zero_reference_metrics(self):
        rho = np.ones((4, 4, 4))
        a = np.zeros((4, 4, 4, 3, 3))
        tau = mixed_stress(rho, rho, 1, a)
        expected = np.broadcast_to(-(2 / 3) * np.eye(3), tau.shape)
        np.testing.assert_allclose(tau, expected, atol=TOL, rtol=0)
        np.testing.assert_array_equal(mixed_stress(rho, np.zeros_like(rho), 1, a), 0)
        metrics = compare_closure(rho, np.zeros(rho.shape + (3,)), 3, 0.25)
        self.assertIsNone(metrics["stress"]["relative_rms_error"])
        self.assertEqual(metrics["stress"]["rms_error"], 0)

    def test_validation_and_honest_synthetic_label(self):
        values = np.ones((4, 4, 4))
        for width in (0, 5, 2.5, True):
            with self.assertRaises(ValueError):
                box_filter(values, width)
        with self.assertRaises(ValueError):
            favre_reference(-values, np.zeros((4, 4, 4, 3)), 3, 0.25)
        with self.assertRaises(ValueError):
            mixed_stress(values, -values, 1, np.zeros((4, 4, 4, 3, 3)))
        result = synthetic_metrics(8, 3)
        self.assertFalse(result["validated_turbulent_reference"])
        self.assertLess(result["trace_identity_max_error"], TOL)
        self.assertIn("mean_signed_error", result["production"])

    def test_signed_and_density_metrics_hand_case(self):
        exact = np.array([-2., -1., 1., 2.])
        density = np.array([1., 2., 3., 4.])
        metrics = production_metrics(exact, np.array([-1., -2., 2., 1.]), density)
        self.assertAlmostEqual(metrics["signed_correlation"], 0.8, delta=TOL)
        self.assertEqual(metrics["exact_backscatter_fraction"], 0.5)
        self.assertEqual(metrics["false_sign_fraction_nonzero_exact"], 0)
        metrics = production_metrics(exact, np.array([1., -2., 2., -1.]), density)
        self.assertAlmostEqual(metrics["signed_correlation"], 0, delta=TOL)
        self.assertEqual(metrics["false_sign_fraction_nonzero_exact"], 0.5)
        low = metrics["subsets"]["density_at_or_below_median"]
        high = metrics["subsets"]["density_above_median"]
        self.assertEqual(low["sample_count"], 2)
        self.assertEqual(high["sample_count"], 2)
        self.assertAlmostEqual(low["errors"]["mean_signed_error"], 1, delta=TOL)
        self.assertAlmostEqual(high["errors"]["mean_signed_error"], -1, delta=TOL)
        self.assertAlmostEqual(low["errors"]["rms_error"], np.sqrt(5), delta=TOL)
        zeros = production_metrics(np.zeros(2), np.ones(2), np.ones(2))
        self.assertIsNone(zeros["signed_correlation"])
        self.assertIsNone(zeros["false_sign_fraction_nonzero_exact"])
        self.assertIsNone(zeros["subsets"]["density_above_median"]["errors"])


if __name__ == "__main__":
    unittest.main()
