"""Exact rational checks for the output-aware selector run envelope.

Fix rational lambda>1 and0<x<=1. Solve for normalized state marker s so
M=[[1,s*x],[x,s*(x+(1-x)/(q-1))]] has Perron eigenvalue lambda.
The actual state marker is y=s/(q-1). This parameterization permits
exact run and envelope checks without numerical eigenvalues.

The helper does not evaluate a global geometry moment or search markers.
It checks the C-retaining envelope and its volume-only limitation.
Only the parent runs --checks; no work is launched on import.
"""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from operator import index
import os
import shutil
import threading
import time


_OS_GF4 = ((0, 0, 0, 0), (0, 1, 2, 3),
           (0, 2, 3, 1), (0, 3, 1, 2))


def output_selector_parameters(q, x, lam):
    """Exact Perron data, with rational lambda chosen as an input."""
    q, x, lam = index(q), Fraction(x), Fraction(lam)
    if q < 2 or not 0 < x <= 1 or lam <= 1:
        raise ValueError('Require q>=2,0<x<=1,lambda>1')
    diagonal = x + (1 - x) / (q - 1)
    s = lam * (lam - 1) / ((lam - 1) * diagonal + x * x)
    y = s / (q - 1)
    rho = (lam - 1) / (s * x)
    kappa = 1 + s * diagonal - lam
    d = (1 + s * x * x) / lam
    c = max(Fraction(1), 1 / rho)
    matrix = ((Fraction(1), s * x), (x, s * diagonal))
    assert 1 + s * x * rho == lam
    assert x + s * diagonal * rho == lam * rho
    assert 0 <= kappa < 1 < lam and x <= rho and 0 < d <= 1
    assert lam / s >= x and x / rho == 1 - kappa
    alpha = lam * (1 - kappa) / (lam - kappa)
    assert d - alpha == kappa * kappa * (lam - 1) / (lam * (lam - kappa))
    assert d >= alpha
    return dict(q=q, x=x, s=s, y=y, lam=lam, rho=rho,
                kappa=kappa, d=d, c=c, matrix=matrix)


def output_selector_run_transfer(parameters, length, forced, closed):
    """Exact bivariate run value; a forced start has a fixed nonzero value."""
    length = index(length)
    if length < 1:
        raise ValueError('Positive run length required')
    matrix = parameters['matrix']
    row = (Fraction(0), Fraction(1)) if forced else (Fraction(1), Fraction(0))
    for _ in range(length - int(forced)):
        row = (row[0] * matrix[0][0] + row[1] * matrix[1][0],
               row[0] * matrix[0][1] + row[1] * matrix[1][1])
    value = row[0] + row[1] * (parameters['x'] if closed else 1)
    return value * (parameters['x'] * parameters['y'] if forced else 1)


def _os_mul(a, b, q):
    if q == 2:
        return a * b
    if q == 4:
        return _OS_GF4[a][b]
    raise ValueError('Amplitude enumeration uses only GF(2) and GF(4)')


@lru_cache(maxsize=256)
def _os_run_histogram(q, delta, length, forced, closed, birth):
    """Independent histogram of(actual state volume,actual output weight)."""
    histogram = Counter()
    for tail in product(range(q), repeat=length - int(forced)):
        values = (birth,) + tail if forced else tail
        previous, weight, volume = 0, 0, 0
        for value in values:
            volume += value != 0
            weight += (value ^ _os_mul(delta, previous, q)) != 0
            previous = value
        if closed:
            weight += _os_mul(delta, previous, q) != 0
        histogram[(volume, weight)] += 1
    return tuple((volume, weight, count) for (volume, weight), count in histogram.items())


def output_selector_run_enumeration(parameters, delta, length, forced, closed, birth=1):
    q, x, y = parameters['q'], parameters['x'], parameters['y']
    if q not in (2, 4) or length < 1 or not 0 < delta < q or (forced and not 0 < birth < q):
        raise ValueError('Invalid small-field run parameters')
    return sum(count * y ** volume * x ** weight
               for volume, weight, count in _os_run_histogram(q, delta, length, forced, closed, birth))


def output_selector_envelope_checks():
    run_checks = volume_controls = death_hole_checks = production_checks = 0
    cases = ((2, (3, 1, 2, 2), (1, 0, 2, 0)),
             (2, (1, 2, 3, 1, 1), (1, 2, 0, 0, 0)),
             (3, (5, 3, 6, 2), (5, 2, 0, 0)))
    for q in (2, 4):
        for x, lam in product((Fraction(1, 4), Fraction(1, 2), Fraction(1)),
                              (Fraction(5, 4), Fraction(2), Fraction(4))):
            data = output_selector_parameters(q, x, lam)
            y, rho, d, c, s = (data[key] for key in ('y', 'rho', 'd', 'c', 's'))
            for length in range(1, 5):
                for forced, closed in product((False, True), repeat=2):
                    transfer = output_selector_run_transfer(data, length, forced, closed)
                    bound = lam ** (length - int(forced))
                    if forced:
                        bound *= x * y * rho
                    bound *= d if closed else c
                    assert transfer <= bound
                    for delta in range(1, q):
                        for birth in (range(1, q) if forced else (1,)):
                            assert output_selector_run_enumeration(data, delta, length, forced, closed, birth) == transfer
                            run_checks += 1
            for volume, births, W, width in ((8, 2, 12, 4), (10, 1, 12, 4), (12, 1, 13, 4)):
                m = volume - births
                ratio = x ** (-W) * (lam / s) ** m * (x * rho) ** births * c ** width
                assert ratio >= x ** (volume + births - W) >= 1
                volume_controls += 1
            for width, allowed, births in cases:
                optional = [(t, i) for t, mask in enumerate(allowed)
                            for i in range(width) if (mask & ~births[t]) >> i & 1]
                f = sum(mask.bit_count() for mask in births)
                allowed_volume = sum(mask.bit_count() for mask in allowed)
                A0 = allowed_volume - f
                closed = sum((allowed[t] & ~allowed[t + 1]).bit_count()
                             for t in range(len(allowed) - 1))
                for bits in product((0, 1), repeat=len(optional)):
                    support = list(births)
                    for value, (t, i) in zip(bits, optional):
                        support[t] |= value << i
                    volume = sum(mask.bit_count() for mask in support)
                    deaths = sum((support[t] & ~support[t + 1]).bit_count()
                                 for t in range(len(support) - 1))
                    holes, m = allowed_volume - volume, volume - f
                    assert deaths <= closed + holes
                    assert lam ** A0 * d ** closed <= (lam / d) ** A0 * d ** (deaths + m)
                    death_hole_checks += 1
    for x, lam in ((Fraction(1, 256), Fraction(2)),
                   (Fraction(1, 65536), Fraction(65536)),
                   (Fraction(1, 1 << 32), Fraction(2)),
                   (Fraction(1), Fraction(4))):
        output_selector_parameters(1 << 32, x, lam)
        production_checks += 1
    print(dict(exact_bivariate_run_checks=run_checks,
               volume_only_negative_controls=volume_controls,
               death_hole_envelope_checks=death_hole_checks,
               production_field_algebra_checks=production_checks,
               global_geometry_sum_bounded=False, full_code_certificate=False), flush=True)


def output_selector_envelope_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true', required=True)
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        output_selector_envelope_checks()
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    output_selector_envelope_main()
