"""GF16/t16 proof development at additional source and state sizes.

This module generalizes the preserved K16 local proof to an explicitly
recorded state dimension. The state6 profile is a different construction,
not an unchanged transfer of the K16 certificate. Floating searches below
are not certificates. No implementation timing is performed here.
"""

import argparse
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
import json
import math

import numpy as np

from gf16_certificate import (CONFIG, HIGH, LOW, Interval,
                              upper_polynomial_product, upper_exact_source, upper_coupled)
from batched_transfer import multiply
from q1_batched_screen import fixed_support_multiplicity
from state_screen import birth_zero_count, zero_patterns, local as generic_local
from subfield_certificate import integer_margin, ZERO


@lru_cache(maxsize=None)
def singleton_weights(state):
    if not 1 <= state <= 8:
        raise ValueError('CA=0 is used only for monomial state sizes through eight')
    weights = []
    for j in range(16):
        weight = 0
        for lane in range(16):
            product = multiply(j, lane, 0x13)
            power, value = 1, int(j == lane)
            for _ in range(state):
                value ^= power
                power = multiply(power, product, 0x13)
            weight += value != 0
        weights.append(weight)
    return tuple(weights)


def local_model(z, state=4):
    result = generic_local(z, 16, state, 16)
    denominator = 16**state-1
    joint = sum(z**w for w in singleton_weights(state))/(16*denominator)
    result[1, 1, 1] += z**16/denominator-joint
    result[1, 1, 0] = joint
    for w in range(2, 17):
        inverse = sum(math.comb(16, h)*z**(16-h)*(1-z)**h
                      /15**max(0, h-16+w) for h in range(17))/denominator
        result[w, 1, 0] = min(result[w, 1, 0], inverse)
    return result


def interval_local_envelope(z, state=4):
    """Same monomial A/C construction, with the explicitly recorded state size."""
    weights = singleton_weights(state)
    if z.lo <= 0 or z.hi > 1:
        raise ValueError('require 0 < z <= 1')
    f = Interval.rational(1, 15)+Interval.rational(14, 15)*z
    zs = [z**i for i in range(17)]
    fs = [f**i for i in range(17)]
    complements = [(1-z)**i for i in range(17)]
    patterns = zero_patterns(16, state, 16)
    denominator = 16**state-1
    result = []
    for w in range(17):
        p0 = Interval.rational(birth_zero_count(w, state, 16), 15**w)
        g = ZERO
        for zeros, count in enumerate(patterns):
            for active_zeros in range(max(0, zeros-(16-w)), min(w, zeros)+1):
                inactive_zeros = zeros-active_zeros
                coefficient = count*math.comb(w, active_zeros)*math.comb(16-w, inactive_zeros)
                g += (Interval.rational(coefficient, denominator)
                      *zs[16-w-inactive_zeros+active_zeros]*fs[w-active_zeros])
        ret = (1-p0)/denominator
        jlo = ret*zs[16]
        jhi = ret*zs[max(0, 17-state-w)]
        if w == 1:
            jlo = jhi = sum((zs[v] for v in weights), ZERO)/(16*denominator)
        elif w >= 2:
            inverse = sum((Interval.rational(math.comb(16, h), denominator*15**max(0, h-16+w))
                           *zs[16-h]*complements[h] for h in range(17)), ZERO)
            jhi = Interval(min(jhi.lo, inverse.lo), min(jhi.hi, inverse.hi))
        result.append([[zs[w]*p0, zs[w]*(1-p0)], [jhi, g-jlo]])
    return result


def geometry(k_log):
    if k_log not in (16, 18, 20):
        raise ValueError('this checkpoint supports K16, K18, and K20 only')
    k = 1 << k_log
    rows = k//8
    n = 16*rows
    return rows, n, (n+9)//10-1


def source_eta(rho):
    """Exact positive affine envelope for every possible nonzero outer support."""
    rho = Fraction(rho)
    if rho <= 0:
        raise ValueError('the source tilt must be positive')
    return max(Fraction(fixed_support_multiplicity(16, 8, w, 16), 16**(w-8))/rho**w
               for w in range(9, 17))


def interval_fraction(value):
    value = Fraction(value)
    return Interval.rational(value.numerator, value.denominator)


def upper_tilted_source(regions, active_rows, rho):
    rho = Fraction(rho)
    eta = source_eta(rho)
    rhoi = interval_fraction(rho)
    matrix = [Decimal(0)]*4
    for h, regional in enumerate(regions):
        weight = HIGH.multiply(Decimal(math.comb(active_rows, h)*16**h), (rhoi**h).hi)
        matrix = [HIGH.add(a, HIGH.multiply(weight, b)) for a, b in zip(matrix, regional)]
    vector = [Decimal(1), Decimal(0)]
    for _ in range(16):
        vector = [HIGH.add(HIGH.multiply(vector[0], matrix[0]), HIGH.multiply(vector[1], matrix[2])),
                  HIGH.add(HIGH.multiply(vector[0], matrix[1]), HIGH.multiply(vector[1], matrix[3]))]
    total = HIGH.multiply(HIGH.add(*vector), (interval_fraction(eta)**active_rows).hi)
    return HIGH.divide(total, Decimal(16**(8*active_rows)))


def upper_regions(z, maximum_h, rows, state=4):
    if rows % 16 or not 1 <= maximum_h < rows:
        raise ValueError('invalid regional geometry')
    local = interval_local_envelope(z, state)
    factors = [[HIGH.multiply(entry.hi, Decimal(math.comb(16, w)))
                for row in matrix for entry in row]
               for w, matrix in enumerate(local[:maximum_h+1])]
    out = [[Decimal(1), Decimal(0), Decimal(0), Decimal(1)]]
    power = rows//16
    while power:
        if power & 1:
            out = upper_polynomial_product(out, factors, maximum_h)
        power >>= 1
        if power:
            factors = upper_polynomial_product(factors, factors, maximum_h)
    return [[HIGH.divide(entry, Decimal(math.comb(rows, h))) for entry in matrix]
            for h, matrix in enumerate(out)]


def class_upper(k_log, active_rows, numerator, denominator, exact_source=None, state=4, source_rho=None):
    rows, n, threshold = geometry(k_log)
    if exact_source is None:
        exact_source = active_rows <= 4
    z = Interval.rational(numerator, denominator)
    regions = upper_regions(z, active_rows, rows, state)
    moment = (upper_exact_source(regions, active_rows) if exact_source
              else upper_tilted_source(regions, active_rows, source_rho) if source_rho is not None
              else upper_coupled(regions, active_rows))
    top = HIGH.multiply(moment, Decimal(math.comb(rows, active_rows)))
    bottom = LOW.multiply(Decimal(15), (z**threshold).lo)
    bound = HIGH.divide(top, bottom)
    return dict(k_log=k_log, state=state, active_rows=active_rows, z=[numerator, denominator],
                source_rho=str(Fraction(source_rho)) if source_rho is not None else None,
                method='exact-source' if exact_source else 'coupled-source-envelope',
                failure_upper=str(bound), integer_bits=integer_margin(bound))


def screen_bound(k_log, active_rows, z, exact_source=None, state=4, source_rho=None):
    from reduced_distance_certificate import log_polynomial_product
    from scipy.special import logsumexp
    rows, n, threshold = geometry(k_log)
    if exact_source is None:
        exact_source = active_rows <= 4
    with np.errstate(divide='ignore'):
        local = np.log(local_model(z, state)[:active_rows+1])
    for w in range(len(local)):
        local[w] += math.log(math.comb(16, w))
    regions = np.array([[[0, -math.inf], [-math.inf, 0]]])
    power = rows//16
    while power:
        if power & 1:
            regions = log_polynomial_product(regions, local, active_rows)
        power >>= 1
        if power:
            local = log_polynomial_product(local, local, active_rows)
    for h in range(active_rows+1):
        regions[h] -= math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
    if exact_source:
        if not 1 <= active_rows <= 4:
            raise ValueError('exact source screen is limited to four rows')
        matrices = [np.exp(regions[h])*16**(h-active_rows/2) for h in range(active_rows+1)]
        masks = [tuple((mask >> j) & 1 for j in range(active_rows)) for mask in range(1 << active_rows)]
        current = np.zeros((1,)*active_rows+(2,))
        current[(0,)*active_rows+(0,)] = 1
        scale = 0.0
        for region in range(16):
            next_state = np.zeros((region+2,)*active_rows+(2,))
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                next_state[target] += current @ matrices[sum(mask)]
            norm = float(np.max(next_state))
            next_state /= norm
            scale += math.log(norm)
            current = next_state
        factors = np.array([fixed_support_multiplicity(16, 8, w, 16)/16**(w-8) if w >= 9 else 0
                            for w in range(17)])
        for axis in range(active_rows):
            shape = [1]*(active_rows+1)
            shape[axis] = 17
            current *= factors.reshape(shape)
        moment_log = scale+math.log(float(np.sum(current)))
    else:
        log_rho = math.log(source_rho) if source_rho is not None else 0.0
        log_eta = (max(math.log(fixed_support_multiplicity(16, 8, w, 16))
                       -(w-8)*math.log(16)-w*log_rho for w in range(9, 17))
                   if source_rho is not None else 0.0)
        source = np.array([math.log(math.comb(active_rows, h))+(h-active_rows/2)*math.log(16)+h*log_rho
                           for h in range(active_rows+1)])
        matrix = logsumexp(regions+source[:, None, None], axis=0)
        vector = np.array([0, -math.inf])
        for region in range(16):
            vector = logsumexp(vector[:, None]+matrix, axis=0)
        moment_log = float(logsumexp(vector))+active_rows*log_eta
    choose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    return -(choose-math.log(15)-threshold*math.log(z)+moment_log)/math.log(2)


def search(k_log, active_rows, exact_source=None, state=4, source_rho=None):
    from scipy.optimize import minimize_scalar
    _, n, _ = geometry(k_log)
    def objective(logtilt):
        z = math.exp(-math.exp(logtilt)*active_rows/n)
        return -screen_bound(k_log, active_rows, z, exact_source, state, source_rho)
    optimum = minimize_scalar(objective, bounds=(1, 8), method='bounded', options=dict(xatol=1e-6))
    return dict(k_log=k_log, state=state, active_rows=active_rows,
                z=math.exp(-math.exp(optimum.x)*active_rows/n), margin_bits=-float(optimum.fun),
                exact_source=(active_rows <= 4 if exact_source is None else exact_source),
                source_rho=source_rho, rigorous=False)


def search_source(k_log, active_rows, state=4):
    from scipy.optimize import minimize
    _, n, _ = geometry(k_log)
    initial = search(k_log, active_rows, exact_source=False, state=state, source_rho=.871)
    logtilt = math.log(-math.log(initial['z'])*n/active_rows)
    def objective(point):
        lt, lr = point
        if not 1 <= lt <= 8 or not -2 <= lr <= 2:
            return 1e100
        z = math.exp(-math.exp(lt)*active_rows/n)
        return -screen_bound(k_log, active_rows, z, False, state, math.exp(lr))
    optimum = minimize(objective, [logtilt, math.log(.871)], method='Nelder-Mead',
                       options=dict(maxiter=160, xatol=1e-6, fatol=1e-6))
    return dict(k_log=k_log, state=state, active_rows=active_rows,
                z=math.exp(-math.exp(optimum.x[0])*active_rows/n), margin_bits=-float(optimum.fun),
                source_rho=math.exp(optimum.x[1]), rigorous=False)


def middle_screen(k_log, state, active, x, z, source_rho=Fraction(871, 1000)):
    from coupled_chernoff import log_terminal_power
    rows, n, threshold = geometry(k_log)
    rho, eta = float(source_rho), float(source_eta(source_rho))
    local = local_model(z, state)
    matrix = sum(math.comb(16, w)*(16*rho*x)**w*(1+x)**(16-w)*local[w] for w in range(17))
    choose = math.lgamma(rows+1)-math.lgamma(active+1)-math.lgamma(rows-active+1)
    value = (-15*choose-8*active*math.log(16)+active*math.log(eta)-math.log(15)
             -threshold*math.log(z)-16*active*math.log(x)+log_terminal_power(matrix, rows))
    return -value/math.log(2)


def middle_search(k_log, state, active, source_rho=Fraction(871, 1000)):
    from scipy.optimize import minimize
    rows, n, threshold = geometry(k_log)
    def objective(point):
        logx, logz = point
        if not -35 < logx < 10 or not -20 < logz < 0:
            return 1e100
        return -middle_screen(k_log, state, active, math.exp(logx), math.exp(logz), source_rho)
    best = None
    for initial_z in [.5, .9, .99, .999]:
        initial_x = active/(max(1, rows-active)*16)
        optimum = minimize(objective, [math.log(initial_x), math.log(initial_z)], method='Nelder-Mead',
                           options=dict(maxiter=300, xatol=1e-9, fatol=1e-7))
        candidate = dict(k_log=k_log, state=state, active_rows=active,
                         x=math.exp(optimum.x[0]), z=math.exp(optimum.x[1]),
                         source_rho=str(source_rho), margin_bits=-float(optimum.fun), rigorous=False)
        if best is None or candidate['margin_bits'] > best['margin_bits']:
            best = candidate
    return best


def middle_upper(k_log, state, active, x, z, source_rho):
    from coupled_chernoff import interval_bound
    rows, n, threshold = geometry(k_log)
    rhoi = interval_fraction(source_rho)
    def local_provider(zi):
        return [[[entry*rhoi**w for entry in row] for row in matrix]
                for w, matrix in enumerate(interval_local_envelope(zi, state))]
    return (interval_bound(active, x, z, local_provider, rows=rows)
            *interval_fraction(source_eta(source_rho))**active)


def propose(k_log=18, state=6, middle_first=224, dense_first=23992):
    """Find rational witnesses; this output has no certificate claim."""
    sparse = []
    margin_sum = 0.0
    for active in range(1, middle_first):
        if 5 <= active <= 16:
            best = search_source(k_log, active, state)
            rho = Fraction(round(best['source_rho']*10**6), 10**6)
        else:
            rho = None if active <= 4 else Fraction(871, 1000)
            best = search(k_log, active, state=state, source_rho=float(rho) if rho else None)
        z = Fraction(round(best['z']*2**40), 2**40)
        bits = screen_bound(k_log, active, float(z), state=state, source_rho=float(rho) if rho else None)
        margin_sum += 2**(-bits)
        sparse.append(dict(active_rows=active, z=str(z), source_rho=str(rho) if rho else None))
    middle = []
    first = middle_first
    while first < dense_first:
        best = middle_search(k_log, state, first)
        x = Fraction(round(best['x']*10**12), 10**12)
        z = Fraction(round(best['z']*10**12), 10**12)
        rho = Fraction(871, 1000)
        if middle_screen(k_log, state, first, float(x), float(z), rho) < 80:
            raise ValueError(f'insufficient middle margin at a={first}')
        low, high = first, dense_first-1
        while low < high:
            probe = (low+high+1)//2
            if middle_screen(k_log, state, probe, float(x), float(z), rho) >= 80:
                low = probe
            else:
                high = probe-1
        middle.append(dict(first=first, last=low, x=str(x), z=str(z), source_rho=str(rho)))
        first = low+1
    return dict(k_log=k_log, state=state, sparse=sparse, middle=middle,
                dense_first=dense_first, dense_z='1/52', dense_bits=80,
                sparse_screen_bits=-math.log2(margin_sum), rigorous=False)


def configuration(k_log, state):
    rows, n, threshold = geometry(k_log)
    singleton_weights(state)
    return dict(CONFIG, k=1 << k_log, n=n, state=state,
                bad_weight_threshold=threshold, refresh=f'uniform-gl{state}-pre-output',
                route='independent-uniform-per-region', initial_state='zero',
                final_state='discarded', regions=16, rows_per_region=rows)


def validate_witnesses(witnesses):
    """Require an exact partition of every possible nonzero outer-row count."""
    rows, _, _ = geometry(witnesses['k_log'])
    singleton_weights(witnesses['state'])
    expected = 1
    if not witnesses['sparse'] or not witnesses['middle']:
        raise ValueError('sparse and middle coverage must both be nonempty')
    for active, z, rho in witnesses['sparse']:
        if active != expected or not 0 < Fraction(z) < 1:
            raise ValueError('invalid sparse coverage or tilt')
        if (active <= 4) != (rho is None):
            raise ValueError('exact source counts are required through four rows')
        if rho is not None and Fraction(rho) <= 0:
            raise ValueError('invalid source tilt')
        expected += 1
    for first, last, x, z, rho in witnesses['middle']:
        if (first != expected or last < first or Fraction(x) <= 0
                or not 0 < Fraction(z) < 1 or Fraction(rho) <= 0):
            raise ValueError('invalid middle coverage or tilt')
        expected = last+1
    if (expected != witnesses['dense_first'] or not expected <= rows
            or not 0 < Fraction(witnesses['dense_z']) < 1
            or not isinstance(witnesses['dense_bits'], int) or witnesses['dense_bits'] < 0):
        raise ValueError('invalid dense-tail coverage or witness')
    return rows


def certify(witnesses):
    """Replay rational witnesses using directed arithmetic, without searches."""
    from subfield_bound_audit import exact_bulk_tail_check
    rows = validate_witnesses(witnesses)
    k_log, state = witnesses['k_log'], witnesses['state']
    sparse, sparse_total = [], Decimal(0)
    for active, ztext, rho in witnesses['sparse']:
        z = Fraction(ztext)
        entry = class_upper(k_log, active, z.numerator, z.denominator,
                            state=state, source_rho=rho)
        sparse.append(entry)
        sparse_total = HIGH.add(sparse_total, Decimal(entry['failure_upper']))
    middle, middle_total = [], Decimal(0)
    for first, last, xtext, ztext, rho in witnesses['middle']:
        x, z = Fraction(xtext), Fraction(ztext)
        left = middle_upper(k_log, state, first, x, z, rho).hi
        right = middle_upper(k_log, state, last, x, z, rho).hi
        contribution = HIGH.multiply(Decimal(last-first+1), max(left, right))
        middle.append(dict(first=first, last=last, x=xtext, z=ztext, source_rho=rho,
                           endpoint_upper=[str(left), str(right)],
                           contribution_upper=str(contribution)))
        middle_total = HIGH.add(middle_total, contribution)
    dense = exact_bulk_tail_check(first=witnesses['dense_first'],
                                  z=Fraction(witnesses['dense_z']), q=16,
                                  k=1 << k_log, margin_bits=witnesses['dense_bits'])
    if not dense['rigorous']:
        raise ValueError('the exact dense-tail comparison failed')
    dense_upper = Interval.rational(1, 1 << witnesses['dense_bits']).hi
    total = HIGH.add(HIGH.add(sparse_total, middle_total), dense_upper)
    return dict(format='gf16-sized-full-v1', complete=True,
                config=configuration(k_log, state), covered_active_rows=[1, rows],
                witnesses=witnesses, sparse=dict(entries=sparse, failure_upper=str(sparse_total)),
                middle=dict(entries=middle, failure_upper=str(middle_total)), dense=dense,
                failure_upper=str(total), integer_bits=integer_margin(total))


def verify(record):
    """Regenerate every claimed contribution and reject incomplete coverage."""
    witnesses = record['witnesses']
    rows = validate_witnesses(witnesses)
    if (record.get('format') != 'gf16-sized-full-v1' or record.get('complete') is not True
            or record.get('config') != configuration(witnesses['k_log'], witnesses['state'])
            or record.get('covered_active_rows') != [1, rows]):
        raise ValueError('incorrect full certificate scope')
    actual = certify(witnesses)
    exact_total = Fraction(0)
    for name, bound_key in [('sparse', 'failure_upper'), ('middle', 'contribution_upper')]:
        supplied_entries, actual_entries = record[name]['entries'], actual[name]['entries']
        if len(supplied_entries) != len(actual_entries):
            raise ValueError('missing or duplicated contribution')
        partial = Fraction(0)
        for supplied, replay in zip(supplied_entries, actual_entries):
            for key, value in replay.items():
                if key not in (bound_key, 'integer_bits', 'endpoint_upper') and supplied.get(key) != value:
                    raise ValueError('incorrect contribution parameters')
            upper = Decimal(supplied[bound_key])
            if not upper.is_finite() or upper < Decimal(replay[bound_key]):
                raise ValueError('understated contribution')
            if 'integer_bits' in supplied and supplied['integer_bits'] > integer_margin(upper):
                raise ValueError('overstated contribution margin')
            if 'endpoint_upper' in replay:
                verify_middle_entry(supplied, replay, upper)
            partial += Fraction(upper)
        upper = Decimal(record[name]['failure_upper'])
        if not upper.is_finite() or Fraction(upper) < partial:
            raise ValueError('understated block sum')
        exact_total += Fraction(upper)
    if record['dense'] != actual['dense']:
        raise ValueError('incorrect exact dense-tail comparison')
    exact_total += Fraction(1, 1 << witnesses['dense_bits'])
    upper = Decimal(record['failure_upper'])
    if not upper.is_finite() or Fraction(upper) < exact_total:
        raise ValueError('understated complete sum')
    if record['integer_bits'] > integer_margin(upper):
        raise ValueError('overstated complete margin')
    return True


def verify_middle_entry(supplied, replay, contribution):
    endpoints = supplied.get('endpoint_upper', [])
    if len(endpoints) != 2:
        raise ValueError('missing middle endpoint')
    values = []
    for given, expected in zip(endpoints, replay['endpoint_upper']):
        value = Decimal(given)
        if not value.is_finite() or value < Decimal(expected):
            raise ValueError('understated middle endpoint')
        values.append(Fraction(value))
    if Fraction(contribution) < (replay['last']-replay['first']+1)*max(values):
        raise ValueError('middle contribution does not cover its displayed endpoints')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, choices=[16, 18, 20], default=18)
    parser.add_argument('--rows', nargs='+', type=int, default=[1, 2, 3, 4, 8, 16, 24, 32, 64, 128])
    parser.add_argument('--state', type=int, choices=list(range(1, 9)), default=4)
    parser.add_argument('--source-rho', type=float)
    parser.add_argument('--optimize-source', action='store_true')
    parser.add_argument('--middle', action='store_true')
    parser.add_argument('--propose', action='store_true')
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--verify')
    args = parser.parse_args()
    if args.verify:
        with open(args.verify, encoding='utf-8-sig') as source:
            verify(json.load(source))
        print('Verified all active-row classes and every rational witness.')
        raise SystemExit(0)
    if args.certify:
        from gf16_size_witnesses import K18_S6
        print(json.dumps(certify(K18_S6), indent=2))
        raise SystemExit(0)
    if args.propose:
        print(json.dumps(propose(args.k_log, args.state), indent=2))
        raise SystemExit(0)
    for active in args.rows:
        answer = (middle_search(args.k_log, args.state, active) if args.middle
                  else search_source(args.k_log, active, args.state) if args.optimize_source
                  else search(args.k_log, active, state=args.state, source_rho=args.source_rho))
        print(json.dumps(answer), flush=True)
