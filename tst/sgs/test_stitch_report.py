"""Independent discontinuity and known time-series fixtures for restart reports."""
import unittest
import numpy as np
from stitch_report import INVENTORIES, join, stationarity
from paper_plots import summarize


class StitchTests(unittest.TestCase):
    def row(self, t):
        return dict(time=t, **{key: 2*t+10 for key in INVENTORIES}, extra=9)

    def test_boundary_once_and_counter_reset_rejected(self):
        origin = self.row(2)
        rows = join([self.row(1), self.row(2)], [self.row(2), self.row(3)], origin)
        self.assertEqual([r['time'] for r in rows], [1, 2, 3])
        self.assertEqual(rows[1]['extra'], 9)
        for key in INVENTORIES:
            reset = dict(origin, **{key: 0})
            with self.assertRaises(ValueError):
                join([self.row(1), self.row(2)], [self.row(3)], reset)

    def test_missing_boundary_and_overlapping_time_rejected(self):
        with self.assertRaises(ValueError):
            join([self.row(1)], [self.row(3)], self.row(2))
        with self.assertRaises(ValueError):
            join([self.row(1), self.row(2)], [self.row(1.5)], self.row(2))

    def data(self, drift):
        ts = np.array([0., .1, 1., 2.3, 5., 7., 10.])
        samples = [dict(time=t, M_rms=5+drift*t, M_sgs=1., K=.1,
                        rho_epsilon=.5, epsilon=.5, rho_mean=1., cs2_mean=1.) for t in ts]
        hist = [dict(time=t, **{'tot-E': 100+2*t, 'represented-E': 100+2*t,
                               'force-work': 2*t, 'sgs-production-integral': t,
                               'sgs-dissipation-integral': t}) for t in ts]
        return dict(protocol=dict(T=1., L=1., V=1., n=64), samples=samples,
                    history=hist, numerical_status='pass')

    def test_stationary_constant_and_known_drift(self):
        protocol = dict(quantities=['M_rms', 'M_sgs', 'K', 'rho_epsilon'],
                        maximum_absolute_fractional_linear_drift=.1,
                        maximum_early_late_relative_difference=.1,
                        maximum_production_dissipation_imbalance=.05)
        for drift, expected in [(0., 'pass'), (.25, 'fail')]:
            data = self.data(drift)
            result = summarize(data)
            self.assertTrue(result['complete_paper_window'])
            self.assertAlmostEqual(result['statistics']['M_rms']['mean'], 5+6*drift)
            screen = stationarity(data, result, protocol)
            self.assertEqual(screen['status'], expected)
            self.assertAlmostEqual(
                screen['metrics']['M_rms']['early_late_relative_difference'],
                4*abs(drift)/(5+6*drift))
        data = self.data(0.)
        for row in data['history']:
            row['sgs-production-integral'] *= 1.052
        self.assertEqual(stationarity(data, summarize(data), protocol)['status'], 'fail')
        data['samples'] = data['samples'][:-1]
        self.assertEqual(
            stationarity(data, summarize(data), protocol)['status'], 'incomplete')


if __name__ == '__main__':
    unittest.main()
