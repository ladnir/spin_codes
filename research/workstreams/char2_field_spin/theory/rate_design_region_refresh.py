"""One-row screen for GL refresh once per region, not once per inner step.

Each region uses a fixed GF16 monomial width-16 A0,C0=A0^T. At region
entry, independently sample a uniform GL(s,256) basis for the global
state. The local region executes Y=X+A0*S, S'=S+C0*X, with no refresh
inside that region. Source symbols have independent nonzero GF256 input
multipliers. The outer is fixed GRS[b,k] over GF256, and there is one
regional global routing stage.

This file bounds only ONE active outer row. It cannot certify the code.
The nonreturn moment uses the deterministic MDS output-distance bound;
all floating optimization is exploratory. No native timings are run.
"""

import argparse
from fractions import Fraction
import json
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from large_state_screen import mul
from q1_batched_screen import fixed_support_multiplicity


def geometric_mean(log_value, count):
    if log_value == 0:
        return 1.
    return -math.expm1(count*log_value)/(-math.expm1(log_value)*count)


def singleton_weights(state):
    emitted, returned = [], []
    for j in range(16):
        values = []
        for lane in range(16):
            product = mul(j, lane)
            power, value = 1, 0
            for _ in range(state):
                value ^= power
                power = mul(power, product)
            values.append(value)
        emitted.append(sum(value != 0 for value in values))
        returned.append(sum((value ^ int(lane == j)) != 0 for lane, value in enumerate(values)))
    return emitted, returned


def region_matrices(z, steps, state, q=256):
    if not 0 < z <= 1 or not 1 <= state <= 8 or steps < 1:
        raise ValueError('invalid one-row regional geometry')
    minimum, denominator = 17-state, q**state-1
    logz = math.log(z)
    empty_nonzero = sum(math.comb(16, w)*fixed_support_multiplicity(16, state, w, q)
                        /denominator*math.exp(logz*steps*w) for w in range(minimum, 17))
    emitted, returned = singleton_weights(state)
    birth = sum(z*geometric_mean(logz*w, steps) for w in emitted)/16
    joint = sum(z**v*geometric_mean(logz*w, steps) for w, v in zip(emitted, returned))/(16*denominator)
    continued = (1-1/denominator)*z**(steps*minimum-1)
    return np.array([[1., 0.], [0., empty_nonzero]]), np.array([[0., birth], [joint, continued]])


def one_row_screen(k_log, z, length=32, dimension=24, state=3):
    multiple = dimension*16
    k = (((1 << k_log)+multiple-1)//multiple)*multiple
    rows, n = k//dimension, k*length//dimension
    threshold = (n+9)//10-1
    r0, r1 = region_matrices(z, rows//16, state)
    with np.errstate(divide='ignore'):
        r0, r1 = np.log(r0), np.log(r1)
    values = np.array([[0., -math.inf]])
    for _ in range(length):
        next_values = np.full((len(values)+1, 2), -math.inf)
        next_values[:-1] = logsumexp(values[:, :, None]+r0, axis=1)
        next_values[1:] = np.logaddexp(next_values[1:], logsumexp(values[:, :, None]+r1, axis=1))
        values = next_values
    contributions = []
    for weight in range(length-dimension+1, length+1):
        count = fixed_support_multiplicity(length, dimension, weight, 256)
        contributions.append(math.log(rows)+math.log(count)-math.log(255)
                             +float(logsumexp(values[weight]))-threshold*math.log(z))
    return -float(logsumexp(contributions))/math.log(2)


def search(k_log=20, length=32, dimension=24, state=3):
    multiple = dimension*16
    k = (((1 << k_log)+multiple-1)//multiple)*multiple
    n = k*length//dimension
    def objective(tilt):
        return -one_row_screen(k_log, math.exp(-math.exp(tilt)/n), length, dimension, state)
    result = minimize_scalar(objective, bounds=(-2, 14), method='bounded')
    return dict(k=k, n=n, outer=[length, dimension], state=state, q=256,
                refresh='independent-uniform-GL-once-per-region', scope='one-active-row-only',
                z=math.exp(-math.exp(result.x)/n), margin_bits=-float(result.fun), rigorous=False)


def return_count_upper(k_log=20, length=32, dimension=24, state=4, q=256, period=16):
    """Exact rational union bound for all one-row words of the new ensemble."""
    if period % 16:
        raise ValueError('injection period must be a multiple of inner width')
    multiple = dimension*period
    k = (((1 << k_log)+multiple-1)//multiple)*multiple
    rows, n = k//dimension, k*length//dimension
    threshold = (n+9)//10-1
    per_region = (17-state)*(rows//16-1)
    permitted = threshold//per_region
    denominator, total, entries = q**state-1, Fraction(0), []
    for weight in range(length-dimension+1, length+1):
        # Every return needs an earlier birth. Return positions cannot be
        # consecutive in active-region order or occupy the first position.
        returns = max(0, (weight-permitted)//2)
        patterns = math.comb(weight-returns, returns)
        words = math.comb(length, weight)*fixed_support_multiplicity(length, dimension, weight, q)//(q-1)
        contribution = Fraction(rows*words*patterns, denominator**returns)
        entries.append(dict(weight=weight, required_returns=returns, return_patterns=patterns,
                            projective_words=words, contribution=contribution))
        total += contribution
    return dict(k=k, n=n, rows=rows, outer=[length, dimension], state=state, q=q, period=period,
                scope='all-one-active-row-words-only', rigorous=True, per_intermediate_region=per_region,
                bad_weight_threshold=threshold, permitted_intermediate_regions=permitted,
                failure_upper=total, entries=entries)


def window_firstmoment_obstruction(k_log=20, period=16, active=12, length=32, dimension=24, state=4, q=256):
    """Exact lower bound on expected bad projective words, not failure probability."""
    multiple = dimension*period
    k = (((1 << k_log)+multiple-1)//multiple)*multiple
    rows, n = k//dimension, k*length//dimension
    threshold = (n+9)//10-1
    window = threshold//(length*16)*16
    kernel_dimension = dimension*active-state*(length-1)
    if kernel_dimension <= 0 or not active <= window < rows:
        raise ValueError('invalid positive-dimensional window witness')
    probability = Fraction(math.comb(window, active), math.comb(rows, active))**length
    lines = (q**kernel_dimension-1)//(q-1)
    return dict(k=k, n=n, rows=rows, active=active, window=window, total_window=window*length,
                bad_weight_threshold=threshold, kernel_dimension=kernel_dimension,
                route_probability=probability, projective_words_lower=lines,
                expected_bad_words_lower=probability*lines, is_failure_probability=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--length', type=int, default=32)
    parser.add_argument('--dimension', type=int, default=24)
    parser.add_argument('--state', type=int, nargs='+', default=[3, 4])
    args = parser.parse_args()
    for state in args.state:
        print(json.dumps(search(args.k_log, args.length, args.dimension, state)), flush=True)
