"""Exact all-h budgets for adaptive state/output zeros in live quartets.

Inject quartet_selector_capacity.py definitions for its directed logs.
For a fixed route, partitions and outer message supported on h rows,
every fresh live quartet gives four independent uniform q-symbol states.
For fixed prescribed coordinate positions, each required to be live and
zero, the joint cost is at most q^-k; requiring each to be live and to
have either state or output zero costs at most(2/q)^k.  Consequently

  Pr[E>=e] <= gamma*C(M,h)*q^(128h)*C(N,e)*q^-e,
  Pr[K_bad>=k] <= gamma*C(M,h)*q^(128h)*C(N,k)*(2/q)^k.

Here E=F-Ls exactly, and K_bad counts each coordinate once.  All later
input activations remain in the full process; extinction is not absorbing.
The H/message union is taken once, with no support-history union.

The coarse thresholds256h and272h each give a global class below2^-239.
Tighter per-h thresholds use certified binary-entropy upper bounds and
allocate2^-92 perh.  Bisection is restricted to1..N/2; beyond that range
the constant binomial bound2^N gives a separate safe linear threshold.
Every returned threshold is checked directly with exact scaled integers.
The output is a partial-class certificate, not a full distance proof.
Only the parent runs jobs; no files, subprocesses or numerical witnesses.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time


AZ_N = 1 << 21
AZ_M = 8192
AZ_Q = 1 << 32
AZ_W = 943718
AZ_SCALE = 1 << 48
AZ_TAIL_BITS = 92
AZ_HVALUES = (1, 4, 16, 64, 128, 256, 512, 768, 1024, 1536, 2048,
              2560, 3072, 3584, 4096, 5120, 6144, 7168)


def _az_ceildiv(a, b):
    if b <= 0:
        raise ValueError('Positive divisor required')
    return -((-a) // b)


def _az_dependencies():
    if globals().get('SC_SCALE') != AZ_SCALE or globals().get('SC_TERMS') != 20:
        raise RuntimeError('Inject the fixed2^48/20-term selector log helpers')
    for name in ('selector_log2_interval', 'selector_binomial_log_interval'):
        if name not in globals():
            raise RuntimeError('Inject quartet_selector_capacity.py directed log helpers')
    # For every mantissa in[1,2], the atanh tail after20 terms is<1
    # scaled unit.  Floor losses are at most20 units.  With ln2>1/2,
    # the normalized log interval has width<128 scaled units.
    if not 2 * AZ_SCALE * 9 < 41 * 8 * 3 ** 41:
        raise ArithmeticError('Uniform atanh-tail comparison failed')
    ln2_lower, ln2_upper = globals()['_SC_LN2']
    if not AZ_SCALE // 2 < ln2_lower <= ln2_upper < AZ_SCALE:
        raise ArithmeticError('The ln2 denominator interval is invalid')
    if ln2_upper - ln2_lower > 21 or 1024 * AZ_N >= AZ_SCALE:
        raise ArithmeticError('Uniform log-width margin failed')


def adaptive_zero_entropy_interval(n, k):
    """Scaled directed interval for n*H2(k/n), not a floating log-binomial."""
    n, k = index(n), index(k)
    if n < 1 or not 0 <= k <= n:
        raise ValueError('Require n>=1 and0<=k<=n')
    if k in (0, n):
        return 0, 0
    lp = selector_log2_interval(k, n)
    lq = selector_log2_interval(n - k, n)
    if lp[1] - lp[0] >= 1024 or lq[1] - lq[0] >= 1024:
        raise ArithmeticError('Directed log interval is unexpectedly wide')
    lower = -k * lp[1] - (n - k) * lq[1]
    upper = -k * lp[0] - (n - k) * lq[0]
    if not 0 <= upper - lower < 1024 * n:
        raise ArithmeticError('Entropy interval width failed')
    return lower, upper


def adaptive_zero_threshold(h, kind='state', tail_bits=AZ_TAIL_BITS):
    """Exact fixed-message class; state rate32, joint-bad-coordinate rate31."""
    _az_dependencies()
    h, tail_bits = index(h), index(tail_bits)
    if not 1 <= h <= AZ_M or tail_bits < 1 or kind not in ('state', 'joint'):
        raise ValueError('Require1<=h<=8192, positive budget and state/joint kind')
    rate = 32 if kind == 'state' else 31
    choose_upper = selector_binomial_log_interval(h)[1]
    # gamma<2 is charged once; q^(128h) counts all outer messages onH.
    base = (4096 * h + 1) * AZ_SCALE + choose_upper
    target = -tail_bits * AZ_SCALE

    def upper_at(e):
        return base + adaptive_zero_entropy_interval(AZ_N, e)[1] - rate * e * AZ_SCALE

    half = AZ_N // 2
    if upper_at(half) <= target:
        low, high = 1, half
        while low < high:
            middle = (low + high) // 2
            if upper_at(middle) <= target:
                high = middle
            else:
                low = middle + 1
        threshold = low
        exponent_upper = upper_at(threshold)
        previous = upper_at(threshold - 1) if threshold > 1 else None
        if previous is not None and previous <= target:
            raise ArithmeticError('Entropy-threshold predecessor still qualifies')
        method = 'directed_entropy_low_half'
        least_for_this_entropy_upper = True
    else:
        threshold = max(half + 1,
                        _az_ceildiv(base + (AZ_N + tail_bits) * AZ_SCALE,
                                    rate * AZ_SCALE))
        exponent_upper = base + AZ_N * AZ_SCALE - rate * threshold * AZ_SCALE
        previous = None
        method = 'constant_binomial_high_half'
        least_for_this_entropy_upper = False
    if threshold > AZ_N:
        return dict(h=h, kind=kind, threshold=AZ_N + 1,
                    method='deterministically_impossible', tail_bits=tail_bits,
                    exact_integer_certificate=True, full_code_certificate=False)
    if exponent_upper > target:
        raise ArithmeticError('Direct threshold certification failed')
    return dict(h=h, kind=kind, threshold=threshold, rate=rate,
                method=method, tail_bits=tail_bits,
                exact_log2_tail_upper=_az_ceildiv(exponent_upper, AZ_SCALE),
                log2_binomial_H_upper=choose_upper,
                least_for_reported_entropy_upper=least_for_this_entropy_upper,
                directed_log_scale=AZ_SCALE,
                exact_integer_certificate=True, no_state_history_union=True,
                future_input_restarts_retained=True, accepted_density_charged_once=True,
                full_code_certificate=False)


def adaptive_zero_threshold_table(hvalues=AZ_HVALUES, tail_bits=AZ_TAIL_BITS):
    return {(index(h), kind): adaptive_zero_threshold(h, kind, tail_bits)
            for h in hvalues for kind in ('state', 'joint')}


def adaptive_zero_profile(h, R, Ls, factual, F, death_lower=0, table=None):
    """Necessary bad-output zero counts; an optional supplied death lower bound."""
    h, R, Ls, factual, F, death_lower = map(index, (h, R, Ls, factual, F, death_lower))
    if not 1 <= h <= AZ_M or not 0 <= factual <= Ls <= F <= AZ_N or death_lower < 0:
        raise ValueError('Invalid actual-state profile')
    state = adaptive_zero_threshold(h, 'state') if table is None else table[(h, 'state')]
    joint = adaptive_zero_threshold(h, 'joint') if table is None else table[(h, 'joint')]
    E = F - Ls
    joint_lower = max(E, F + death_lower - AZ_W)
    return dict(h=h, R=R, Ls=Ls, factual=factual, F=F, E=E,
                death_lower=death_lower, required_joint_bad_coordinates=joint_lower,
                state_zero_threshold=state['threshold'], joint_bad_threshold=joint['threshold'],
                state_zero_class=E >= state['threshold'],
                joint_bad_class=joint_lower >= joint['threshold'],
                complement_death_upper=AZ_W - F + joint['threshold'] - 1,
                full_code_certificate=False)


def adaptive_zero_coarse_constants():
    """Integer/rational constants in the two uniform all-h<2^-239 classes."""
    assert AZ_N == 256 * AZ_M
    assert 3 * AZ_N < (1 << 15) * 256
    assert 3 * AZ_N < (1 << 15) * 272
    assert 3 * AZ_M < 1 << 15
    assert 32 * 128 - (32 - 15) * 256 + 15 == -241
    assert 32 * 128 - (31 - 15) * 272 + 15 == -241
    # Sum_h>=1 2^(1-241h) <2^-239.
    assert 2 ** 241 > 2
    assert Fraction(2, (1 << 241) - 1) < Fraction(1, 1 << 239)
    assert AZ_M * 2 == 1 << 14
    assert 14 - AZ_TAIL_BITS == -78
    return dict(state_coarse_multiple=256, joint_coarse_multiple=272,
                each_coarse_global_exponent=-239,
                each_tight_all_h_exponent=-79,
                both_tight_all_h_exponent=-78,
                full_code_certificate=False)


def adaptive_zero_checks():
    """Exact small binomial and overlap identities; parent-only execution."""
    _az_dependencies()
    constants = adaptive_zero_coarse_constants()
    entropy_checks = zero_union_checks = overlap_checks = 0
    for n in range(2, 65):
        previous = None
        for k in range(1, n + 1):
            _, upper = adaptive_zero_entropy_interval(n, k)
            assert selector_log2_interval(comb(n, k))[1] <= upper + 1024
            if k <= n // 2:
                value = upper - 31 * k * AZ_SCALE
                if previous is not None:
                    assert value < previous
                previous = value
            entropy_checks += 1
    # A fixed-coordinate bad event has one forbidden q-value when old=0,
    # and two when old!=0.  Count a whole quartet without independence claims.
    for q in (4, 8):
        for old in ((0, 0, 0, 0), (1, 0, 1, 0), (1, 1, 1, 1)):
            state = joint = shared = 0
            from itertools import product
            for new in product(range(q), repeat=4):
                qzero = [value == 0 for value in new]
                yzero = [value == before for value, before in zip(new, old)]
                state += all(qzero)
                joint += all(a or b for a, b in zip(qzero, yzero))
                shared += sum(a and b for a, b in zip(qzero, yzero))
                E = sum(qzero)
                survivor_zeros = sum(y and before != 0 for y, before in zip(yzero, old))
                assert sum(a or b for a, b in zip(qzero, yzero)) == E + survivor_zeros
                overlap_checks += 1
            assert state == 1 and joint <= 16
            zero_union_checks += 1
    for h in (1, 4, 256, 4096, 8192):
        for kind in ('state', 'joint'):
            row = adaptive_zero_threshold(h, kind)
            assert row['threshold'] <= AZ_N
            assert row['exact_log2_tail_upper'] <= -AZ_TAIL_BITS
    print(dict(exact_entropy_checks=entropy_checks,
               exact_quartet_zero_union_checks=zero_union_checks,
               exact_state_output_overlap_checks=overlap_checks,
               threshold_rows=10, **constants), flush=True)


def adaptive_zero_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--table', action='store_true')
    parser.add_argument('--hvalues', type=int, nargs='+', default=AZ_HVALUES)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            adaptive_zero_checks()
        else:
            print(adaptive_zero_coarse_constants(), flush=True)
            for row in adaptive_zero_threshold_table(args.hvalues).values():
                print(row, flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_zero_main()
