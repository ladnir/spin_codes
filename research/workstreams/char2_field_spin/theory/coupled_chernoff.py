"""Coefficient-Chernoff bounds for the one-route regional field construction.

The optimization is a floating screen, not a certificate. ``interval_bound``
accepts an independently outward-rounded local transition provider and
rational witnesses; it does not use floating arithmetic.
"""

import argparse
from fractions import Fraction
import json
import math

import numpy as np

from subfield_certificate import Interval, ONE, ZERO, interval_matrix_product


# Search witnesses only: interval_cover rechecks their mathematical claims.
GF16_K16_MIDDLE_WITNESSES = (
    dict(first=128, last=191, x='1204610467/200000000000', z='905721316261/1000000000000'),
    dict(first=192, last=336, x='9233110209/1000000000000', z='862077669351/1000000000000'),
    dict(first=337, last=670, x='3322935107/200000000000', z='96818722129/125000000000'),
    dict(first=671, last=1238, x='30067254791/1000000000000', z='642977317237/1000000000000'),
    dict(first=1239, last=1650, x='8786946783/200000000000', z='64329513083/125000000000'),
    dict(first=1651, last=2333, x='38238979947/500000000000', z='131898833471/500000000000'),
    dict(first=2334, last=5963, x='224852056129/1000000000000', z='78963199943/1000000000000'),
    dict(first=5964, last=5997, x='449428897973/250000000000', z='14030444437/1000000000000'),
)


def log_terminal_power(matrix, exponent):
    """Log of e_zero * matrix**exponent * (1,1), with positive scaling."""
    matrix_scale = float(np.max(matrix))
    matrix = matrix / matrix_scale
    log_scale = math.log(matrix_scale)
    vector = np.array([1.0, 0.0])
    accumulated = 0.0
    while exponent:
        if exponent & 1:
            vector = vector @ matrix
            scale = float(np.max(vector))
            vector /= scale
            accumulated += log_scale + math.log(scale)
        exponent >>= 1
        if exponent:
            matrix = matrix @ matrix
            scale = float(np.max(matrix))
            matrix /= scale
            log_scale = 2*log_scale + math.log(scale)
    return accumulated + math.log(float(np.sum(vector)))


def screen_bound(active, x, z, rows=8192, q=16, state=4,
                 width=16, length=16, dimension=8, distance=0.1):
    from reduced_distance_certificate import local_model
    if not 1 <= active <= rows or rows % width:
        raise ValueError('invalid regional geometry')
    local = local_model(z, width, q, state)
    matrix = sum(math.comb(width, w)*(q*x)**w*(1+x)**(width-w)*local[w]
                 for w in range(width+1))
    logchoose = (math.lgamma(rows+1)-math.lgamma(active+1)
                 - math.lgamma(rows-active+1))
    threshold = math.ceil(distance*rows*length)-1
    value = (-(length-1)*logchoose-(length-dimension)*active*math.log(q)
             - math.log(q-1)-threshold*math.log(z)-length*active*math.log(x)
             + log_terminal_power(matrix, rows*length//width))
    return -value/math.log(2)


def search(active, **parameters):
    from scipy.optimize import minimize
    best = None

    def objective(point):
        logx, logz = point
        if not (-35 < logx < 10 and -20 < logz < 0):
            return 1e100
        return -screen_bound(active, math.exp(logx), math.exp(logz), **parameters)

    rows, q = parameters.get('rows', 8192), parameters.get('q', 16)
    for initial_z in [0.3, 0.7, 0.9, 0.98, 0.999]:
        initial_x = active/(max(1, rows-active)*q)
        result = minimize(objective, [math.log(initial_x), math.log(initial_z)],
                          method='Nelder-Mead',
                          options=dict(maxiter=400, xatol=1e-10, fatol=1e-8))
        if best is None or -result.fun > best['margin_bits']:
            best = dict(active_rows=active, x=math.exp(result.x[0]),
                        z=math.exp(result.x[1]), margin_bits=-float(result.fun),
                        rigorous=False, method='coupled-coefficient-chernoff')
    return best


def interval_bound(active, x, z, local_provider, rows=8192, q=16,
                   width=16, length=16, dimension=8, distance=Fraction(1, 10)):
    """Outward bound for rational positive x and 0<z<1.

    ``local_provider`` must return entrywise moment upper bounds, enclosed
    in Interval entries, for all weights zero through width at this z.
    The bound covers this entire active-row class, including outer-word
    multiplicities and one projective division by q-1.
    """
    x, z = Fraction(x), Fraction(z)
    if not 1 <= active <= rows or rows % width or x <= 0 or not 0 < z < 1:
        raise ValueError('invalid certificate parameters')
    xi = Interval.rational(x.numerator, x.denominator)
    zi = Interval.rational(z.numerator, z.denominator)
    local = local_provider(zi)
    matrix = [[ZERO, ZERO], [ZERO, ZERO]]
    for w in range(width+1):
        factor = math.comb(width, w)*(q*xi)**w*(ONE+xi)**(width-w)
        for i in range(2):
            for j in range(2):
                matrix[i][j] += factor*local[w][i][j]
    out = [[ONE, ZERO], [ZERO, ONE]]
    exponent = rows*length//width
    while exponent:
        if exponent & 1:
            out = interval_matrix_product(out, matrix)
        exponent >>= 1
        if exponent:
            matrix = interval_matrix_product(matrix, matrix)
    desired = Fraction(distance)*rows*length
    threshold = (desired.numerator+desired.denominator-1)//desired.denominator-1
    denominator = (Interval.integer(math.comb(rows, active))**(length-1)
                   * Interval.integer(q)**((length-dimension)*active)
                   * (q-1)*zi**threshold*xi**(length*active))
    return (out[0][0]+out[0][1])/denominator


def search_cover(first=128, last=5997, class_margin=80):
    """Choose rational witnesses; every resulting range still needs checking.

    For fixed x,z the logarithm of the class bound is convex in a, because
    binom(L,a) is log-concave and every other a-dependent factor is geometric.
    Therefore checking both endpoints bounds all classes in the interval.
    """
    cover = []
    while first <= last:
        best = search(first)
        x = Fraction(round(best['x']*10**12), 10**12)
        z = Fraction(round(best['z']*10**12), 10**12)
        if screen_bound(first, float(x), float(z)) < class_margin:
            raise ValueError(f'no adequate witness at a={first}')
        low, high = first, last
        while low < high:
            middle = (low+high+1)//2
            if screen_bound(middle, float(x), float(z)) >= class_margin:
                low = middle
            else:
                high = middle-1
        cover.append(dict(first=first, last=low, x=str(x), z=str(z)))
        first = low+1
    return cover


def interval_cover(witnesses, local_provider, first=128, last=5997):
    """Outward sum over every class in consecutive witnessed intervals."""
    expected, total = first, ZERO
    details = []
    for witness in witnesses:
        begin, end = witness['first'], witness['last']
        if begin != expected or end < begin or end > last:
            raise ValueError('the cover has a gap, overlap, or invalid endpoint')
        endpoints = [interval_bound(a, witness['x'], witness['z'], local_provider)
                     for a in (begin, end)]
        upper = max(bound.hi for bound in endpoints)
        contribution = Interval(upper, upper)*(end-begin+1)
        total += contribution
        details.append(dict(**witness, contribution_upper=str(contribution.hi)))
        expected = end+1
    if expected != last+1:
        raise ValueError('the cover does not reach the requested final class')
    return total, details


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', nargs='+', type=int, default=[64, 128, 256, 512, 1024, 4096, 5997])
    parser.add_argument('--cover', action='store_true')
    parser.add_argument('--verify-middle', action='store_true')
    args = parser.parse_args()
    if args.verify_middle:
        from gf16_certificate import interval_local_envelope
        from subfield_certificate import integer_margin
        total, details = interval_cover(GF16_K16_MIDDLE_WITNESSES, interval_local_envelope)
        print(json.dumps(dict(covered_active_rows=[128, 5997], failure_upper=str(total.hi),
                              integer_bits=integer_margin(total.hi), rigorous=True,
                              complete=False, witnesses=details), indent=2))
    elif args.cover:
        print(json.dumps(search_cover(), indent=2), flush=True)
    else:
        for active in args.active:
            print(json.dumps(search(active)), flush=True)
