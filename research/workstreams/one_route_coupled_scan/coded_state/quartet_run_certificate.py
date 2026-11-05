"""Exact positive-vector certificate for the run/occupancy moment.

The proposed vector is floating, but every row of K v <= lambda v is
verified with integers.  The frozen parameters are s=75/64, z^4=1/2,
t=d(1/256)=(q+65536)/(256(q+257)), q=2^32.  Lambda is a power of two.
The free-endpoint bound is lambda^T*v[0]/min(v), with no terminal reset.

The exact verification uses two positive triangular integer transforms;
it constructs no dense rational matrix and no large field powers.  No
witness files or dumps are written.  This certifies only the geometry
moment, not the whole-code distance theorem.  Run all jobs serially.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb
import os
import shutil
import threading
import time

import numpy as np


RC_LD = np.longdouble
RC_Q = 1 << 32
RC_S = Fraction(75, 64)
RC_T = Fraction(RC_Q + 65536, 256 * (RC_Q + 257))


def run_certificate_ceil_log2_ratio(numerator, denominator):
    if numerator <= 0 or denominator <= 0:
        raise ValueError('Positive integer ratio required')
    exponent = numerator.bit_length() - denominator.bit_length()
    if exponent >= 0:
        return exponent + int(numerator > (denominator << exponent))
    return exponent + int((numerator << -exponent) > denominator)


def run_certificate_occupancy(width):
    """Integer numerators/denominators of M_u, with z^4=1/2."""
    if not isinstance(width, int) or width < 4 or width % 4:
        raise ValueError('Width must be a positive integer multiple of four')
    groups = width // 4
    numerator = [0] * (width + 1)
    polynomial = [1]
    for k in range(groups + 1):
        multiplier = comb(groups, k) << (groups - k)
        for u, coefficient in enumerate(polynomial):
            numerator[u] += multiplier * coefficient
        if k < groups:
            after = [0] * (len(polynomial) + 4)
            for u, coefficient in enumerate(polynomial):
                for j in range(1, 5):
                    after[u + j] += coefficient * comb(4, j)
            polynomial = after
    denominator = [comb(width, u) << groups for u in range(width + 1)]
    return numerator, denominator


def run_certificate_exact_action(vector, s=RC_S, t=RC_T, occupancy=None):
    """Return exact integer numerator/denominator pairs for K*vector.

    Write t=a/b and s=c/d.  First form
      N_j=sum_l C(n-j,l)c^l d^(n-j-l) v_(j+l).
    Then T*v at row u has numerator
      sum_j C(u,j)a^(u-j)[c(b-a)]^j N_j
    and denominator d^n b^u.  All coefficients are nonnegative for t<=1.
    """
    if any(not isinstance(value, int) or value <= 0 for value in vector):
        raise ValueError('The witness must consist of positive Python integers')
    width = len(vector) - 1
    s, t = Fraction(s), Fraction(t)
    if s <= 0 or not 0 <= t <= 1:
        raise ValueError('Require s>0 and 0<=t<=1 for the positive transforms')
    a, b, c, d = t.numerator, t.denominator, s.numerator, s.denominator
    if occupancy is None:
        occupancy = run_certificate_occupancy(width)
    m_numerator, m_denominator = occupancy
    if len(m_numerator) != width + 1 or len(m_denominator) != width + 1:
        raise ValueError('Occupancy data has a different width')
    d_powers = [1]
    for _ in range(width):
        d_powers.append(d_powers[-1] * d)
    background = [0] * (width + 1)
    for j in range(width + 1):
        degree = width - j
        coefficient = d_powers[degree]
        total = coefficient * vector[j]
        for ell in range(degree):
            # This division is exact: the result is a binomial monomial.
            coefficient = coefficient * (degree - ell) * c // ((ell + 1) * d)
            total += coefficient * vector[j + ell + 1]
        background[j] = total
    a_power = b_power = 1
    step = c * (b - a)
    numerator, denominator = [], []
    for u in range(width + 1):
        if a == 0:
            transformed = step ** u * background[u]
        else:
            coefficient = a_power
            transformed = coefficient * background[0]
            for j in range(u):
                coefficient = coefficient * (u - j) * step // ((j + 1) * a)
                transformed += coefficient * background[j + 1]
        numerator.append(m_numerator[u] * transformed)
        denominator.append(m_denominator[u] * d_powers[width] * b_power)
        a_power *= a
        b_power *= b
    return numerator, denominator


def run_certificate_verify(vector, updates=4096, s=RC_S, t=RC_T,
                            target_lambda_log=5, max_lambda_log=6):
    """Verify all rows exactly; return only compact certificate metadata."""
    if not isinstance(updates, int) or updates < 0:
        raise ValueError('The update count must be a nonnegative integer')
    numerator, denominator = run_certificate_exact_action(vector, s, t)
    worst = 0
    for u in range(1, len(vector)):
        if numerator[u] * denominator[worst] * vector[worst] > numerator[worst] * denominator[u] * vector[u]:
            worst = u
    required = run_certificate_ceil_log2_ratio(numerator[worst], denominator[worst] * vector[worst])
    chosen = max(target_lambda_log, required)
    if chosen > max_lambda_log:
        return dict(verified=False, required_lambda_log2=required,
                    maximum_allowed_lambda_log2=max_lambda_log, worst_row=worst,
                    full_code_certificate=False)
    # This final full pass is the actual integer Collatz certificate.
    for u in range(len(vector)):
        left, right = numerator[u], denominator[u] * vector[u]
        if chosen >= 0:
            assert left <= (right << chosen)
        else:
            assert (left << -chosen) <= right
    ratio_bits = run_certificate_ceil_log2_ratio(vector[0], min(vector))
    moment_bits = updates * chosen + ratio_bits
    return dict(verified=True, exact_rows=len(vector), width=len(vector) - 1, updates=updates,
                lambda_log2=chosen, target_lambda_log2=target_lambda_log,
                target_achieved=chosen == target_lambda_log, worst_row=worst,
                terminal_ratio_bits=ratio_bits, moment_log2_integer_upper=moment_bits,
                witness_max_bits=max(value.bit_length() for value in vector),
                witness_min_bits=min(vector).bit_length(),
                s_numerator=Fraction(s).numerator, s_denominator=Fraction(s).denominator,
                t_numerator=Fraction(t).numerator, t_denominator=Fraction(t).denominator,
                occupancy_quartet_numerator=1, occupancy_quartet_denominator=2,
                full_code_certificate=False)


def run_certificate_propose(width=512, iterations=256, mantissa_bits=48):
    """Floating power iteration only proposes integers; it proves no bound."""
    if not isinstance(iterations, int) or iterations < 1 or not 8 <= mantissa_bits <= 60:
        raise ValueError('Require positive iterations and 8..60 mantissa bits')
    if 'run_occupancy_matrix' in globals():
        matrix_function = globals()['run_occupancy_matrix']
    else:
        # Ordinary standalone use; concatenated remote loaders inject this.
        from quartet_run_occupancy import run_occupancy_matrix
        matrix_function = run_occupancy_matrix
    s_float = RC_LD(RC_S.numerator) / RC_LD(RC_S.denominator)
    t_float = RC_LD(RC_T.numerator) / RC_LD(RC_T.denominator)
    matrix, scale = matrix_function(width, np.log2(s_float), np.log2(t_float), -RC_LD(1) / 4)
    vector = np.ones(width + 1, dtype=RC_LD)
    for _ in range(iterations):
        vector = matrix @ vector
        maximum = vector.max()
        if maximum <= 0 or not np.isfinite(maximum) or np.any(vector <= 0):
            raise ArithmeticError('Floating proposal overflow/underflow')
        vector /= maximum
    residual = (matrix @ vector) / vector
    float_ratio_log = scale + np.log2(residual.max())
    mantissas, exponents = np.frexp(vector)
    minimum_exponent = int(exponents.min())
    unit = RC_LD(2) ** mantissa_bits
    integers = [int(np.ceil(mantissa * unit)) << (int(exponent) - minimum_exponent)
                for mantissa, exponent in zip(mantissas, exponents)]
    if min(integers) <= 0:
        raise ArithmeticError('Integer witness proposal is not positive')
    return integers, dict(proposal_iterations=iterations, proposal_mantissa_bits=mantissa_bits,
                          float_max_ratio_log2=float(float_ratio_log), floating_proposal_only=True)


def _rc_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << j for j in rest)]]
    raise ValueError('Independent exact checks use width four or eight')


def _rc_reference(width, s, t):
    """Independent physical-mask/partition count, with z^(4k)=2^-k."""
    partitions = _rc_partitions(width)
    result = np.full((width + 1, width + 1), Fraction(0), dtype=object)
    occupancy = []
    for u in range(width + 1):
        old = (1 << u) - 1
        moment = sum(Fraction(1, 1 << sum(bool(quad & old) for quad in partition))
                     for partition in partitions) / len(partitions)
        occupancy.append(moment)
        for new in range(1 << width):
            result[u, new.bit_count()] += moment * s ** new.bit_count() * t ** (old & ~new).bit_count()
    return result, occupancy


def run_certificate_checks():
    F = Fraction
    rows = occupancies = moments = log_checks = 0
    for numerator in range(1, 33):
        for denominator in range(1, 33):
            exponent = run_certificate_ceil_log2_ratio(numerator, denominator)
            scale = F(2) ** exponent
            assert F(numerator, denominator) <= scale
            assert F(numerator, denominator) > scale / 2
            log_checks += 1
    for width in (4, 8):
        for s, t in ((F(5, 4), F(1, 3)), (RC_S, F(3, 4)),
                     (RC_S, RC_T), (F(3, 4), F(0)), (F(3, 4), F(1))):
            reference, exact_occupancy = _rc_reference(width, s, t)
            occupancy = run_certificate_occupancy(width)
            for i in range(width + 1):
                assert F(occupancy[0][i], occupancy[1][i]) == exact_occupancy[i]
                occupancies += 1
            vector = [1 + i * i for i in range(width + 1)]
            numerator, denominator = run_certificate_exact_action(vector, s, t, occupancy)
            expected = reference @ np.array(vector, dtype=object)
            for i in range(width + 1):
                assert F(numerator[i], denominator[i]) == expected[i]
                rows += 1
            certificate = run_certificate_verify(vector, updates=4, s=s, t=t,
                target_lambda_log=0, max_lambda_log=64)
            assert certificate['verified']
            lam = F(2) ** certificate['lambda_log2']
            assert all(expected[i] <= lam * vector[i] for i in range(width + 1))
            row = np.array([F(int(i == 0)) for i in range(width + 1)], dtype=object)
            for updates in range(5):
                assert sum(row) <= lam ** updates * F(vector[0], min(vector))
                moments += 1
                row = row @ reference
    print(dict(exact_integer_action_rows=rows, exact_occupancy_checks=occupancies,
               free_endpoint_inequalities=moments, integer_log_checks=log_checks,
               floating_proposal_used=False, full_code_certificate=False), flush=True)


def run_certificate_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    parser.add_argument('--width', type=int, default=512)
    parser.add_argument('--updates', type=int, default=4096)
    parser.add_argument('--iterations', type=int, default=256)
    parser.add_argument('--mantissa-bits', type=int, default=48)
    parser.add_argument('--lambda-log', type=int, default=5)
    parser.add_argument('--max-lambda-log', type=int, default=6)
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
            run_certificate_checks()
            return
        vector, proposal = run_certificate_propose(args.width, args.iterations, args.mantissa_bits)
        print(dict(stage='integer_witness_proposed', seconds=round(time.monotonic() - started, 3),
                   **proposal), flush=True)
        certificate = run_certificate_verify(vector, args.updates,
            target_lambda_log=args.lambda_log, max_lambda_log=args.max_lambda_log)
        print(dict(seconds=round(time.monotonic() - started, 3), **certificate), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    run_certificate_main()
