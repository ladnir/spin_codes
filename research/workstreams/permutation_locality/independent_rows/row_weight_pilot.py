"""Independent-row Bernoulli conditioning: selected homogeneous weight classes.

Research-only binary64 pilot using the authenticated nine-coordinate cache.
Not an outward certificate, not a cover, and not the eleven-coordinate model.
"""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import comb, log
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bch_joint_support import authenticated_caps
from categorical_model import mix
from categorical_screen import code_hash
from occupancy_memory import TERMINAL
from occupancy_sensitivity import float_placement
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment, log_binomial_mass


def packet_distribution(probabilities):
    values = [1]
    for p in probabilities:
        previous = values
        values = [0]*(len(previous)+1)
        for weight, mass in enumerate(previous):
            values[weight] += mass*(1-p)
            values[weight+1] += mass*p
    return values


def self_test():
    checks = 0
    for probabilities in ((Q(1, 3),)*4, (Q(0), Q(1, 3), Q(1, 2), Q(1)),
                          (Q(1, 5), Q(2, 5), Q(3, 5), Q(4, 5))):
        law = packet_distribution(probabilities)
        direct = [Q(0)]*5
        for mask in product((0, 1), repeat=4):
            mass = Q(1)
            for bit, p in zip(mask, probabilities):
                mass *= p if bit else 1-p
            direct[sum(mask)] += mass
        assert law == direct and sum(law) == 1
        checks += 6
    # Two-column rows: conditioning each row on its weight gives the
    # independent uniform supports used by the row-shuffle construction.
    probabilities = (Q(1, 5), Q(1, 3), Q(1, 2), Q(3, 4))
    for words in product(range(4), repeat=4):
        weights = [word.bit_count() for word in words]
        probability = Q(1)
        event = Q(1)
        target = Q(1)
        for w, p in zip(weights, probabilities):
            row = p**w*(1-p)**(2-w)
            probability *= row
            event *= comb(2, w)*row
            target /= comb(2, w)
        assert probability/event == target
        checks += 1
    print('Exact Poisson-binomial and independent row-conditioning checks:', checks, 'passed', flush=True)


def load_cache(path, tilt):
    with np.load(path, allow_pickle=False) as saved:
        expected = code_hash()
        if str(saved['code_hash']) != expected or str(saved['tilt']) != tilt:
            raise ValueError('cache source hash or tilt does not match')
        data = [(saved[f's{j}'], saved[f't{j}'], saved[f'c{j}']) for j in range(33)]
    assert all(np.isfinite(matrices).all() and (matrices >= 0).all()
               for _, matrices, _ in data)
    print('Authenticated categorical cache source hash:', expected, flush=True)
    return data


def inner_value(data, q, tilt, p, active_rows=4):
    law = np.array(packet_distribution([p]*active_rows+[0.]*(4-active_rows)), dtype=float)
    active = 1-law[0]
    theta = law[1:]/active
    operators = mix(data, theta)
    region = float_placement(operators, q)
    matrix = matrix_for_probabilities(region, [active]*q)
    return log_power_moment(matrix, 256, TERMINAL)+tilt*209715


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, default=64)
    parser.add_argument('--active-rows', type=int, nargs='+', choices=(1, 2, 3, 4), default=[4])
    parser.add_argument('--weights', type=int, nargs='+', default=[48, 64, 80, 96, 112, 128])
    parser.add_argument('--tilt', default='.032')
    parser.add_argument('--cache', type=Path,
                        default=Path('tmp/categorical-q64-pilot-checked.npz'))
    parser.add_argument('--max-evaluations', type=int, default=40)
    parser.add_argument('--self-test-only', action='store_true')
    args = parser.parse_args()
    if not 1 <= args.groups <= 2048 or any(not 1 <= w <= 255 for w in args.weights):
        parser.error('groups must be in [1,2048] and weights in [1,255]')
    self_test()
    if args.self_test_only:
        return
    data = load_cache(args.cache, args.tilt)
    caps = authenticated_caps()
    tilt = float(args.tilt)
    q = args.groups
    memo = {}
    def inner(p, active_rows):
        key = p, active_rows
        if key not in memo:
            memo[key] = inner_value(data, q, tilt, p, active_rows)
        return memo[key]
    for active_rows in args.active_rows:
        for w in args.weights:
            if not caps[w]:
                print('Weight excluded by authenticated spectrum cap:', w, flush=True)
                continue
            outer = log(comb(2048, q))+q*log(comb(4, active_rows))+active_rows*q*log(caps[w])
            def score(p):
                return (inner(p, active_rows)-active_rows*q*log_binomial_mass(256, w, p)+outer)/log(2)
            initial = w/256
            base = score(initial)
            print('BINARY64 ROW-WEIGHT CLASS q/r/w', q, active_rows, w, 'coarse p', initial,
                  'log2 contribution upper', base, flush=True)
            fit = minimize_scalar(lambda x: score(1/(1+np.exp(-x))),
                                  bounds=(-6, 6), method='bounded',
                                  options={'maxiter': args.max_evaluations, 'xatol': 1e-5})
            p = 1/(1+np.exp(-fit.x))
            best = min((base, initial), (fit.fun, p))
            print('BINARY64 OPTIMIZED q/r/w', q, active_rows, w, 'log2 contribution upper', best[0],
                  'p', best[1], 'gain bits', base-best[0], 'converged', fit.success,
                  'evaluations', fit.nfev, flush=True)
    print('Exactly q active groups, each with r nonzero rows of weight w and 4-r zero rows; row positions counted.', flush=True)
    print('Selected classes only. Binary64 nine-coordinate diagnostic, NOT a certificate or full support cover.', flush=True)


if __name__ == '__main__':
    main()
