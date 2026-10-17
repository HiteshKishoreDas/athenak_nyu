"""Independent ideal-gas normal-shock/path fixtures; no detector imports.

PHY01/02/03/05/06: construct forward Rankine-Hugoniot states, not inverse Mach.
An ideal discontinuity crossing is an oracle time, not a detector-mask entry time.
Numerical mask entry/exit brackets must be reported independently.
"""
import json
import math


def planar_rh(mach=2.0, gamma=1.4, rho1=1.0, pressure1=1.0,
              front0=0.3, speed=0.5, downstream_side=-1):
    """Return left/right states for a front x=front0+speed*t.

    downstream_side=-1 means compressed fluid on the left and incoming
    shock-frame flow toward decreasing x. Reflection uses downstream_side=+1.
    """
    if not all(math.isfinite(v) for v in (mach, gamma, rho1, pressure1,
                                         front0, speed)):
        raise ValueError("Finite shock fixture parameters required")
    if mach <= 1 or gamma <= 1 or rho1 <= 0 or pressure1 <= 0:
        raise ValueError("Adiabatic shock requires M>1, gamma>1, rho1,P1>0")
    if downstream_side not in (-1, 1):
        raise ValueError("downstream_side must be -1 or +1")
    sound1 = math.sqrt(gamma * pressure1 / rho1)
    ratio_rho = (gamma + 1) * mach**2 / ((gamma - 1) * mach**2 + 2)
    ratio_pressure = (2 * gamma * mach**2 - gamma + 1) / (gamma + 1)
    relative1 = downstream_side * mach * sound1
    relative2 = relative1 / ratio_rho
    upstream = dict(rho=rho1, pressure=pressure1, velocity=speed + relative1)
    downstream = dict(rho=rho1 * ratio_rho, pressure=pressure1 * ratio_pressure,
                      velocity=speed + relative2)
    return dict(mach=mach, gamma=gamma, front0=front0, speed=speed,
                downstream_side=downstream_side, upstream=upstream,
                downstream=downstream,
                left=downstream if downstream_side == -1 else upstream,
                right=upstream if downstream_side == -1 else downstream,
                density_ratio=ratio_rho, pressure_ratio=ratio_pressure,
                temperature_ratio=ratio_pressure / ratio_rho)


def crossing_time(fixture, particle0):
    """Ideal upstream material crossing; negative value means prehistory."""
    return ((particle0 - fixture['front0']) /
            (fixture['speed'] - fixture['upstream']['velocity']))


def material_position(fixture, particle0, time):
    """Piecewise analytic path through an infinitesimal planar shock."""
    crossing = crossing_time(fixture, particle0)
    if crossing < 0:
        return particle0 + fixture['downstream']['velocity'] * time
    if time <= crossing:
        return particle0 + fixture['upstream']['velocity'] * time
    at_front = fixture['front0'] + fixture['speed'] * crossing
    return at_front + fixture['downstream']['velocity'] * (time - crossing)


def self_check():
    """Check independent conservation laws, reflection and boost invariance."""
    for gamma in (1.4, 5 / 3):
        for mach in (1.01, 2, 10, 100):
            for side in (-1, 1):
                f = planar_rh(mach=mach, gamma=gamma, downstream_side=side)
                a, b = f['upstream'], f['downstream']
                w1, w2 = a['velocity'] - f['speed'], b['velocity'] - f['speed']
                def equal(x, y):
                    assert math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-12)
                equal(a['rho'] * w1, b['rho'] * w2)
                equal(a['pressure'] + a['rho'] * w1**2,
                      b['pressure'] + b['rho'] * w2**2)
                equal(gamma / (gamma - 1) * a['pressure'] / a['rho'] + w1**2 / 2,
                      gamma / (gamma - 1) * b['pressure'] / b['rho'] + w2**2 / 2)
                xp = f['front0'] - side * 0.25
                tc = crossing_time(f, xp)
                equal(material_position(f, xp, tc), f['front0'] + f['speed'] * tc)
                boosted = planar_rh(mach=mach, gamma=gamma,
                                    downstream_side=side, speed=f['speed'] + 7)
                equal(crossing_time(boosted, xp), tc)
                equal(material_position(boosted, xp, tc + 0.01),
                      material_position(f, xp, tc + 0.01) + 7 * (tc + 0.01))
    f = planar_rh()
    assert math.isclose(f['density_ratio'], 8 / 3)
    assert math.isclose(f['pressure_ratio'], 4.5)
    return {'test_ids': ['PHY01', 'PHY02', 'PHY03', 'PHY05', 'PHY06'],
            'oracle_self_check': 'passed', 'fixture': f,
            'particles': [dict(x0=x, exact_ideal_crossing=crossing_time(f, x))
                          for x in (0.45, 0.6, 0.75)]}


if __name__ == '__main__':
    print(json.dumps(self_check(), indent=2))
