"""Independent exact oracles, metrics and their non-production selfchecks.

Run directly with --output /absolute/path/results.json. The existing regression
runner can discover run()/analyze(); its build is unnecessary for these checks.
All positions are continuous lifts, not minimum-image periodic coordinates.
These prescribed velocities are kinematic fixtures, not free Euler solutions.
"""

import argparse
import json
import math
from pathlib import Path

try:
    from ..utils import athena
except ImportError:  # Standalone oracle checks need no executable or dependency.
    athena = None


# Frozen before candidate runs. Absolute budgets assume unit characteristic length.
ROUNDING_TOL = 2.0e-12
QUADRATURE_TOL = 1.0e-12
POSITION_BUDGET = 1.0e-4
MIN_RK2_ORDER = 1.8
_results = []


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def translation(x, time, velocity):
    return tuple(a + time*b for a, b in zip(x, velocity))


def acceleration(x, time, velocity, rate):
    return tuple(a + time*b + 0.5*time*time*c
                 for a, b, c in zip(x, velocity, rate))


def sinusoidal(x, time, amplitude, frequency, phase=0.0):
    if frequency == 0.0:
        return translation(x, time, tuple(a*math.sin(phase) for a in amplitude))
    integral = (math.cos(phase) - math.cos(frequency*time + phase))/frequency
    return tuple(a + b*integral for a, b in zip(x, amplitude))


def rotation(x, time, omega, axis=(0.0, 0.0, 1.0), centre=(0.0, 0.0, 0.0)):
    """Rodrigues rotation; also maps tangents by using centre=(0,0,0)."""
    size = norm(axis)
    if size == 0.0:
        raise ValueError('rotation axis is zero')
    axis = tuple(a/size for a in axis)
    v = tuple(a-b for a, b in zip(x, centre))
    cross = (axis[1]*v[2] - axis[2]*v[1],
             axis[2]*v[0] - axis[0]*v[2], axis[0]*v[1] - axis[1]*v[0])
    dot = sum(a*b for a, b in zip(axis, v))
    c, s = math.cos(omega*time), math.sin(omega*time)
    return tuple(o + c*a + s*b + (1-c)*dot*d
                 for o, a, b, d in zip(centre, v, cross, axis))


def strain(x, time, rates):
    return tuple(a*math.exp(time*b) for a, b in zip(x, rates))


def shear(x, time, rate):
    return (x[0] + rate*time*x[1], x[1], x[2])


def reversing_shear(x, time, amplitude, wavenumber, frequency):
    """u=(A sin(k*y) cos(w*t),0,0); return at each 2*pi/w."""
    integral = time if frequency == 0.0 else math.sin(frequency*time)/frequency
    return (x[0] + amplitude*math.sin(wavenumber*x[1])*integral, x[1], x[2])


def reversing_tangent(x, tangent, time, amplitude, wavenumber, frequency):
    integral = time if frequency == 0.0 else math.sin(frequency*time)/frequency
    gradient = amplitude*wavenumber*math.cos(wavenumber*x[1])*integral
    return (tangent[0] + gradient*tangent[1], tangent[1], tangent[2])


def circle(label, radius, centre=(0.0, 0.0, 0.0)):
    return (centre[0] + radius*math.cos(label),
            centre[1] + radius*math.sin(label), centre[2])


def helix(label, radius, pitch):
    return (radius*math.cos(label), radius*math.sin(label), pitch*label)


def circle_geometry(radius, samples):
    if radius <= 0.0 or samples < 3:
        raise ValueError('circle needs positive radius and at least three samples')
    return {'length': 2*math.pi*radius, 'curvature': 1/radius,
            'polygon_length': 2*samples*radius*math.sin(math.pi/samples)}


def helix_geometry(radius, pitch, parameter_extent):
    if radius <= 0.0 or parameter_extent <= 0.0:
        raise ValueError('helix needs positive radius and parameter extent')
    return {'length': parameter_extent*math.hypot(radius, pitch),
            'curvature': radius/(radius*radius + pitch*pitch)}


def polyline_length(points, closed=False, closure_displacement=(0.0, 0.0, 0.0)):
    if len(points) < (3 if closed else 2):
        raise ValueError('insufficient samples')
    ends = list(points[1:])
    if closed:
        ends.append(tuple(a+b for a, b in zip(points[0], closure_displacement)))
    return sum(norm(tuple(b-a for a, b in zip(p, q)))
               for p, q in zip(points, ends))


def curvature(p, q, r):
    """Unsigned inverse circumradius; degenerate triples have no valid value."""
    a = tuple(y-x for x, y in zip(p, q))
    b = tuple(y-x for x, y in zip(q, r))
    c = tuple(y-x for x, y in zip(p, r))
    denominator = norm(a)*norm(b)*norm(c)
    if denominator == 0.0:
        return None
    cross = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
             a[0]*b[1]-a[1]*b[0])
    return 2*norm(cross)/denominator


def local_stretch(tangent, initial_tangent):
    denominator = norm(initial_tangent)
    if denominator == 0.0:
        raise ValueError('zero reference tangent')
    return norm(tangent)/denominator


def integrate(function, lower, upper, tolerance=QUADRATURE_TOL):
    """Adaptive Simpson quadrature, independent of production geometry kernels."""
    def recurse(a, b, fa, fm, fb, whole, tol, depth):
        middle = (a+b)/2
        fl, fr = function((a+middle)/2), function((middle+b)/2)
        left = (middle-a)*(fa + 4*fl + fm)/6
        right = (b-middle)*(fm + 4*fr + fb)/6
        difference = left+right-whole
        if abs(difference) <= 15*tol:
            return left+right+difference/15
        if depth == 0:
            raise ValueError('quadrature did not converge')
        return (recurse(a, middle, fa, fl, fm, left, tol/2, depth-1)
                + recurse(middle, b, fm, fr, fb, right, tol/2, depth-1))
    fa, fm, fb = function(lower), function((lower+upper)/2), function(upper)
    whole = (upper-lower)*(fa+4*fm+fb)/6
    return recurse(lower, upper, fa, fm, fb, whole, tolerance, 24)


def position_error(actual, expected, characteristic_length=1.0):
    if (not actual or len(actual) != len(expected)
            or not math.isfinite(characteristic_length) or characteristic_length <= 0.0):
        raise ValueError('position comparison needs matching nonempty fixed labels')
    errors = []
    for p, q in zip(actual, expected):
        if len(p) != 3 or len(q) != 3:
            raise ValueError('positions must have three components')
        error = norm(tuple(a-b for a, b in zip(p, q)))/characteristic_length
        if not math.isfinite(error):
            raise ValueError('nonfinite position error')
        errors.append(error)
    return max(errors)


def length_error(actual, reference):
    if (not math.isfinite(actual) or not math.isfinite(reference)
            or actual < 0.0 or reference <= 0.0):
        raise ValueError('length comparison needs finite positive reference')
    return abs(actual-reference)/reference


def geometric_error(actual, reference, closed=False, subdivisions=8):
    """Sampled symmetric distance to segments; refine subdivisions independently.

    Both input curves must already be lifted, with the winding closure supplied
    as a repeated displaced endpoint for noncontractible periodic curves.
    """
    if subdivisions < 1 or len(actual) < 2 or len(reference) < 2:
        raise ValueError('geometry comparison needs segments and positive subdivisions')
    if any(len(p) != 3 or any(not math.isfinite(x) for x in p)
           for p in list(actual) + list(reference)):
        raise ValueError('geometry comparison needs finite three-component points')
    def segments(points):
        return list(zip(points, points[1:] + ([points[0]] if closed else [])))
    def distance(p, edge):
        a, b = edge
        v = tuple(y-x for x, y in zip(a, b))
        offset = tuple(y-x for x, y in zip(a, p))
        denominator = sum(x*x for x in v)
        fraction = (sum(x*y for x, y in zip(v, offset))/denominator
                    if denominator > 0 else 0.0)
        fraction = max(0.0, min(1.0, fraction))
        return norm(tuple(y-x-fraction*z for x, y, z in zip(a, p, v)))
    def directed(first, second):
        errors = []
        for a, b in first:
            for index in range(subdivisions+1):
                fraction = index/subdivisions
                point = tuple(x+fraction*(y-x) for x, y in zip(a, b))
                error = min(distance(point, edge) for edge in second)
                if not math.isfinite(error):
                    raise ValueError('nonfinite geometry error')
                errors.append(error)
        return max(errors)
    first, second = segments(list(actual)), segments(list(reference))
    return max(directed(first, second), directed(second, first))


def arc_weights(points, closed=False):
    """Trapezoidal arc-length weights; endpoint half-weights for an open string."""
    weights = [0.0]*len(points)
    if len(points) < 2:
        raise ValueError('arc weights need at least two samples')
    for index in range(len(points) if closed else len(points)-1):
        following = (index+1) % len(points)
        length = norm(tuple(a-b for a, b in zip(points[index], points[following])))
        weights[index] += length/2
        weights[following] += length/2
    return weights


def check_positions(actual, expected, tolerance=POSITION_BUDGET):
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError('invalid position tolerance')
    error = position_error(actual, expected)
    if error > tolerance:
        raise ValueError('position error {} exceeds {}'.format(error, tolerance))
    return error


def check_order(errors, budget=POSITION_BUDGET, minimum=MIN_RK2_ORDER):
    if (not math.isfinite(budget) or budget <= 0.0
            or not math.isfinite(minimum) or minimum <= 0.0):
        raise ValueError('invalid convergence gate')
    if len(errors) < 4 or any(not math.isfinite(e) or e <= 0.0 for e in errors):
        raise ValueError('order gate needs four finite positive errors')
    slopes = [math.log(a/b, 2) for a, b in zip(errors, errors[1:])]
    if any(p < minimum for p in slopes[-2:]) or errors[-1] > budget:
        raise ValueError('convergence/absolute budget failed: {}'.format(slopes))
    return slopes


def rk_reference(velocity, initial, time, steps, order=4):
    """Dense reference/selfcheck only; never called to obtain production answers."""
    x, t, dt = initial, 0.0, time/steps
    def advance(base, derivative, factor):
        return tuple(a + factor*b for a, b in zip(base, derivative))
    for _ in range(steps):
        k1 = velocity(x, t)
        if order == 2:
            k2 = velocity(advance(x, k1, dt), t+dt)
            x = tuple(a+dt*(b+c)/2 for a, b, c in zip(x, k1, k2))
        else:
            k2 = velocity(advance(x, k1, dt/2), t+dt/2)
            k3 = velocity(advance(x, k2, dt/2), t+dt/2)
            k4 = velocity(advance(x, k3, dt), t+dt)
            x = tuple(a+dt*(b+2*c+2*d+e)/6
                      for a, b, c, d, e in zip(x, k1, k2, k3, k4))
        t += dt
    return x


def selfchecks():
    records = []
    def check(case, name, actual, expected, tolerance=ROUNDING_TOL):
        error = abs(actual-expected)
        if not math.isfinite(error) or error > tolerance:
            raise ValueError('{} failed: {}'.format(name, error))
        records.append({'case': case, 'name': name, 'status': 'pass',
                        'error': error, 'tolerance': tolerance})
    def rejects(name, function):
        try:
            function()
        except ValueError:
            records.append({'case': 'R04', 'name': name, 'status': 'pass'})
        else:
            raise ValueError('wrong output was accepted: '+name)

    initial = (0.4, -0.3, 0.2)
    check('A02', 'translation', translation(initial, 2, (1, -2, 3))[1], -4.3)
    check('A03', 'acceleration', acceleration(initial, 2, (1, 0, 0),
                                             (2, 0, 0))[0], 6.4)
    check('A03', 'sinusoidal', sinusoidal(initial, math.pi, (1, 0, 0), 1)[0], 2.4)
    check('A04', 'quarter_rotation', rotation((1, 0, 0), math.pi/2, 1)[1], 1)
    oblique = rotation(initial, 0.7, 1.3, (1, 2, 3))
    check('A04', 'oblique_rotation_preserves_norm', norm(oblique), norm(initial))
    check('A05', 'diagonal_strain', strain((1, 2, 0), 1, (1, -1, 0))[0], math.e)
    check('A06', 'shear', shear(initial, 2, 3)[0], -1.4)
    check('A05', 'material_stretch', local_stretch(strain((1, 0, 0), 1,
                                                       (1, -1, 0)), (1, 0, 0)),
          math.e)
    for samples in (8, 16, 32, 64):
        points = [circle(2*math.pi*i/samples, 0.7) for i in range(samples)]
        exact = circle_geometry(0.7, samples)
        check('G02', 'polygon_length_{}'.format(samples),
              polyline_length(points, True), exact['polygon_length'])
        check('G02', 'circle_curvature_{}'.format(samples),
              curvature(*points[:3]), exact['curvature'])
    errors = [length_error(circle_geometry(0.7, n)['polygon_length'], 1.4*math.pi)
              for n in (8, 16, 32, 64)]
    slopes = check_order(errors, budget=5e-4)
    records.append({'case': 'E04', 'name': 'circle_smooth_length_refinement',
                    'status': 'pass', 'errors': errors, 'slopes': slopes,
                    'minimum_order': MIN_RK2_ORDER, 'budget': 5e-4})
    check('G03', 'helix_length', helix_geometry(2, 3, 4)['length'], 4*math.sqrt(13))
    check('G03', 'helix_curvature', helix_geometry(2, 3, 4)['curvature'], 2/13)
    errors, bending_errors = [], []
    for intervals in (8, 16, 32, 64):
        points = [helix(2*math.pi*i/intervals, 0.7, 0.2)
                  for i in range(intervals+1)]
        exact = helix_geometry(0.7, 0.2, 2*math.pi)
        errors.append(length_error(polyline_length(points), exact['length']))
        bending_errors.append(abs(curvature(*points[:3])-exact['curvature']))
    for name, errors, budget in (('helix_length', errors, 5e-4),
                                 ('helix_curvature', bending_errors, 5e-4)):
        slopes = check_order(errors, budget=budget)
        records.append({'case': 'E04', 'name': name+'_refinement', 'status': 'pass',
                        'errors': errors, 'slopes': slopes,
                        'minimum_order': MIN_RK2_ORDER, 'budget': budget})
    check('G01', 'straight_curvature', curvature((0, 0, 0), (1, 0, 0),
                                               (3, 0, 0)), 0)
    check('G04', 'quadrature_circle', integrate(lambda a: 0.7, 0, 2*math.pi),
          1.4*math.pi)
    check('G05', 'nonuniform_arc_weights', sum(arc_weights(
        [(0, 0, 0), (0.1, 0, 0), (3, 0, 0)])), 3)
    check('E04', 'geometric_resampling', geometric_error(
        [(0, 0, 0), (3, 0, 0)], [(0, 0, 0), (0.1, 0, 0), (3, 0, 0)]), 0)
    check('E04', 'geometric_displacement', geometric_error(
        [(0, 0, 0), (3, 0, 0)], [(0, 1, 0), (3, 1, 0)]), 1)
    ellipse = integrate(lambda a: math.hypot(2*math.sin(a), math.cos(a)),
                        0, 2*math.pi)
    refined = integrate(lambda a: math.hypot(2*math.sin(a), math.cos(a)),
                        0, 2*math.pi, QUADRATURE_TOL/16)
    check('E06', 'ellipse_reference_refinement', ellipse, refined, QUADRATURE_TOL)

    time = 0.8
    exact = sinusoidal(initial, time, (1, 0, 0), 1)
    errors = [position_error([rk_reference(lambda x, t: (math.sin(t), 0, 0),
                                          initial, time, n, order=2)], [exact])
              for n in (8, 16, 32, 64)]
    slopes = check_order(errors)
    records.append({'case': 'E02', 'name': 'rk2_reference_harness', 'status': 'pass',
                    'errors': errors, 'slopes': slopes,
                    'minimum_order': MIN_RK2_ORDER, 'budget': POSITION_BUDGET,
                    'coverage': 'reference selfcheck, no production execution'})
    frequency, amplitude, wavenumber = 1.7, 0.8, 2.3
    period = 2*math.pi/frequency
    velocity = lambda x, t: (amplitude*math.sin(wavenumber*x[1])
                             * math.cos(frequency*t), 0, 0)
    exact = reversing_shear(initial, period/4, amplitude, wavenumber, frequency)
    errors = [position_error([rk_reference(velocity, initial, period/4, n)], [exact])
              for n in (4, 8, 16, 32)]
    slopes = check_order(errors, budget=1e-8, minimum=3.8)
    records.append({'case': 'A08', 'name': 'dense_reversing_shear_reference',
                    'status': 'pass', 'errors': errors, 'slopes': slopes,
                    'minimum_order': 3.8, 'budget': 1e-8,
                    'coverage': 'reference convergence, no production execution'})
    check('A08', 'full_period_exact_return', position_error(
        [reversing_shear(initial, period, amplitude, wavenumber, frequency)],
        [initial]), 0)
    check('A08', 'full_period_dense_return', position_error(
        [rk_reference(velocity, initial, period, 128)], [initial]), 0)

    rejects('stale_stage_velocity', lambda: check_positions([initial],
            [acceleration(initial, 1, (0, 0, 0), (1, 0, 0))]))
    rejects('first_order', lambda: check_order([0.008, 0.004, 0.002, 0.001]))
    rejects('one_good_slope', lambda: check_order([1e-4, 8e-5, 7e-5, 1e-5]))
    rejects('nan_output', lambda: check_positions([(math.nan, 0, 0)], [initial]))
    rejects('nan_later_row', lambda: check_positions([initial, (math.nan, 0, 0)],
            [initial, initial]))
    rejects('missing_component', lambda: check_positions([(0, 0)], [(0, 0)]))
    rejects('empty_output', lambda: check_positions([], []))
    rejects('large_absolute_error', lambda: check_order([1, 0.25, 0.0625, 0.015625]))
    rejects('zero_reference_tangent', lambda: local_stretch(initial, (0, 0, 0)))
    points = [circle(2*math.pi*i/16, 0.7) for i in range(16)]
    rejects('missing_closure_length', lambda: check(
        'R04', 'missing_closure', polyline_length(points), polyline_length(points, True)))
    return records


def run(**kwargs):
    global _results
    _results = selfchecks()
    return 'independent oracle selfchecks; executable path not exercised'


def analyze():
    return bool(_results) and all(row['status'] == 'pass' for row in _results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    run()
    result = {'schema': 'athenak-material-curves-oracles-v1',
              'implementation_execution': False, 'status': 'pass' if analyze() else 'fail',
              'thresholds': {'rounding': ROUNDING_TOL, 'quadrature': QUADRATURE_TOL,
                             'position': POSITION_BUDGET, 'rk2_order': MIN_RK2_ORDER},
              'tests': _results}
    report = json.dumps(result, indent=2, allow_nan=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report+'\n')
    else:
        print(report)
    return 0 if analyze() else 1


if __name__ == '__main__':
    raise SystemExit(main())
