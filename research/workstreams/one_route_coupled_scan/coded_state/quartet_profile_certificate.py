"""Exact profile-class certificate from the run/occupancy moment and routing.

Inject quartet_run_occupancy.py and quartet_run_certificate.py definitions
with their main blocks removed.  --certify proposes a vector and REVERIFIES
every row with integers; no previously observed moment bound is hardcoded.
Large binomials, routing exponents, constants and box monotonicity tests use
integers only.  No logarithm or floating comparison contributes to a claim.

The frozen parameters are x=1/256, s=75/64, z^4=1/2 and t=d(x).
For each canonical nonzero profile (h,R,L,f,F), all D are summed by the
moment.  If M<=2^Mbits, a valid integer exponent is

 Mbits-30+8W-32R+24L-floor(117L/512)+(F-f)/4+min(0,routebits).

The routing ratio is C(F-f,C_H-f)C(N-F,256h-C_H)/C(8192,h)^255,
C_H=F-R+128h.  At most 2^96 nonzero profiles exist.  Therefore every
profile satisfying exponent<=-172 belongs to a class whose TOTAL failure
probability is less than 2^-76.  This is a class theorem, not a certificate
for every remaining profile or the full fixed-construction distance goal.

A box uses a uniform PER-PROFILE maximum; the global 96-bit union already
counts every R, so no second R-range union is added.  The free terminal
remains in the freshly verified moment certificate.
"""

import argparse
from fractions import Fraction
from math import comb
import os
import shutil
import threading
import time


PC_N = 1 << 21
PC_K = 1 << 20
PC_WIDTH = 512
PC_T = 4096
PC_ROWS = 8192
PC_REGIONS = 256
PC_W = 943718
PC_Q = 1 << 32
PC_PROFILE_UNION_BITS = 96
PC_THRESHOLD = -172


def profile_certificate_constants():
    q, N = PC_Q, PC_N
    x = Fraction(1, 256)
    c = (q + 1) * x
    lam = 1 + c
    s = Fraction(75, 64)
    t = Fraction(q + 65536, 256 * (q + 257))
    forced_closed = min(Fraction(1), lam / (q + 1))
    unforced_closed = (1 + q * x * x) / lam
    assert t == max(forced_closed, unforced_closed) and 0 < t < 1
    assert c <= lam and lam == (1 << 24) * Fraction(q + 257, q)
    # (1+257/q)^L <= (1+257/q)^N <= 1/(1-257N/q) < 2.
    # The middle inequality follows by bounding each binomial coefficient by N^k.
    assert 257 * N < q
    correction = Fraction(q, q - 257 * N)
    assert correction == Fraction(2048, 1791) and correction < 2
    assert 75 ** 1024 > 1 << 6378
    assert s ** 1024 > 1 << 234  # log2(s)>117/512, without using logarithms.
    assert Fraction(1024, 989) < 2
    assert (1 << 13) * (1 << 20) * ((1 << 21) ** 3) == 1 << PC_PROFILE_UNION_BITS
    assert PC_THRESHOLD + PC_PROFILE_UNION_BITS == -76
    return dict(q=q, x_numerator=1, x_denominator=256,
                s_numerator=75, s_denominator=64,
                closure_numerator=t.numerator, closure_denominator=t.denominator,
                lambda_correction_numerator=correction.numerator,
                lambda_correction_denominator=correction.denominator,
                profile_union_bits=PC_PROFILE_UNION_BITS,
                qualifying_profile_exponent=PC_THRESHOLD,
                resulting_class_failure_exponent=-76)


def profile_certificate_ceil_log2_ratio(numerator, denominator):
    """The smallest integer exponent e with numerator/denominator<=2^e."""
    if numerator <= 0 or denominator <= 0:
        raise ValueError('Positive integer numerator and denominator required')
    exponent = numerator.bit_length() - denominator.bit_length()
    if exponent >= 0:
        return exponent + int(numerator > (denominator << exponent))
    return exponent + int((numerator << -exponent) > denominator)


def _profile_certificate_choose(n, k, cache):
    if n < 0 or k < 0 or k > n:
        raise ValueError('Binomial arguments outside support')
    k = min(k, n - k)
    known = cache.setdefault('binomials', {})
    key = (n, k)
    if key in known:
        return known[key]
    # Reuse nearby point/box binomials by exact short recurrences.  This
    # avoids recomputing a million-bit central binomial for a small box shift.
    choices = []
    for old_n, old_k in known:
        if old_n == n and abs(old_k - k) <= 2048:
            choices.append((abs(old_k - k), old_n, old_k))
        elif old_k == k and abs(old_n - n) <= 2048:
            choices.append((abs(old_n - n), old_n, old_k))
    if choices:
        _, old_n, old_k = min(choices)
        value = known[old_n, old_k]
        if old_n == n:
            while old_k < k:
                value = value * (n - old_k) // (old_k + 1)
                old_k += 1
            while old_k > k:
                value = value * old_k // (n - old_k + 1)
                old_k -= 1
        else:
            while old_n < n:
                value = value * (old_n + 1) // (old_n + 1 - k)
                old_n += 1
            while old_n > n:
                value = value * (old_n - k) // old_n
                old_n -= 1
    else:
        value = comb(n, k)
    known[key] = value
    return value


def profile_certificate_route(h, R, f, F, cache=None):
    """Exact forced-birth routing ratio, including the eligible-set union."""
    if cache is None:
        cache = {}
    C_H = F - R + 128 * h
    selected = 256 * h
    if not (1 <= h <= PC_ROWS and 0 <= f <= C_H <= min(F, selected)
            and 0 <= selected - C_H <= PC_N - F):
        raise ValueError('Routing arguments outside binomial support')
    covered = _profile_certificate_choose(F - f, C_H - f, cache)
    missing = _profile_certificate_choose(PC_N - F, selected - C_H, cache)
    denominators = cache.setdefault('quota_denominators', {})
    if h not in denominators:
        denominators[h] = _profile_certificate_choose(PC_ROWS, h, cache) ** 255
    denominator = denominators[h]
    numerator = covered * missing
    bits = profile_certificate_ceil_log2_ratio(numerator, denominator)
    return dict(routebits=bits, probability_capped_bits=min(0, bits),
                raw_ratio_at_most_one=numerator <= denominator,
                C_H=C_H, eligible_holes=selected - C_H,
                numerator_bit_length=numerator.bit_length(), denominator_bit_length=denominator.bit_length())


def _profile_certificate_exponent(moment_bits, R, L, f, F, routing_bits):
    # For a noncanonical box corner use ceil((F-f)/4); every actual profile
    # has F-f divisible by four, and the continuous monotonicity proof precedes rounding.
    quarter = (F - f + 3) // 4
    return (moment_bits - 30 + 8 * PC_W - 32 * R + 24 * L
            - (117 * L) // 512 + quarter + routing_bits)


def profile_certificate_point(h, R, L, f, F, moment_bits, cache=None):
    if not (1 <= R <= PC_K and 1 <= f <= L <= F <= PC_N):
        raise ValueError('Nonzero profile dimensions are outside their physical ranges')
    if (F - f) % 4:
        raise ValueError('A canonical profile has F-f divisible by four')
    C_H = F - R + 128 * h
    if not 129 * h <= C_H <= 256 * h or R - 128 * h > 128 * (PC_ROWS - h):
        raise ValueError('Profile violates eligible/outside row capacity')
    routing = profile_certificate_route(h, R, f, F, cache)
    exponent = _profile_certificate_exponent(moment_bits, R, L, f, F, routing['probability_capped_bits'])
    return dict(kind='point_profile', h=h, R=R, L=L, f=f, F=F,
                all_availability_closure_counts=True, moment_log2_integer_upper=moment_bits,
                exact_profile_exponent_upper=exponent, qualifies_global_profile_class=exponent <= PC_THRESHOLD,
                qualifying_threshold=PC_THRESHOLD, global_profile_union_bits=PC_PROFILE_UNION_BITS,
                full_code_certificate=False, **routing)


def profile_certificate_box_ratios(N, selected, rank_base, q, R_box, F_box, f_lower):
    """Exact support and adjacent-ratio checks for the coupled weighted numerator."""
    Rlo, Rhi = R_box
    Flo, Fhi = F_box
    Cmin = Flo - Rhi + rank_base
    Cmax = Fhi - Rlo + rank_base
    support = (0 <= Flo <= Fhi <= N and Rlo <= Rhi and 0 <= f_lower <= Cmin
               and Cmax <= selected and 0 <= Rlo - rank_base
               and Rhi - rank_base <= N - selected)
    result = dict(support_valid=support, C_min=Cmin, C_max=Cmax)
    if not support:
        return dict(**result, ratios_valid=False)
    # Largest Q(R+1)/Q(R) in the rectangle occurs at Fhi,Rlo,flo.
    C = Fhi - Rlo + rank_base
    numerator_R = (C - f_lower) * (N - Fhi - selected + C)
    denominator_R = (Fhi - C + 1) * (selected - C + 1)
    R_ok = numerator_R <= q * denominator_R
    # Largest Q(F+1)/Q(F) occurs at Flo,Rhi,flo.  z^4=1/2.
    C = Flo - Rhi + rank_base
    numerator_F = (Flo + 1 - f_lower) * (selected - C)
    denominator_F = (C + 1 - f_lower) * (N - Flo)
    F_ok = Flo == Fhi or (denominator_F > 0 and 2 * numerator_F ** 4 <= denominator_F ** 4)
    return dict(**result, R_ratio_numerator=numerator_R, R_ratio_denominator=denominator_R,
                F_ratio_numerator=numerator_F, F_ratio_denominator=denominator_F,
                R_weighted_ratio_at_most_one=R_ok, F_weighted_ratio_at_most_one=F_ok,
                ratios_valid=R_ok and F_ok)


def profile_certificate_box(h, R_box, L_box, f_box, F_box, moment_bits, cache=None):
    for interval in (R_box, L_box, f_box, F_box):
        if len(interval) != 2 or any(not isinstance(value, int) for value in interval) or interval[0] > interval[1]:
            raise ValueError('Each box coordinate must be a nonempty closed integer interval')
    if not (1 <= h <= PC_ROWS and 1 <= R_box[0] <= R_box[1] <= PC_K
            and all(1 <= lo <= hi <= PC_N for lo, hi in (L_box, f_box, F_box))):
        raise ValueError('Box exceeds the counted nonzero-profile ranges')
    checks = profile_certificate_box_ratios(PC_N, 256 * h, 128 * h, PC_Q, R_box, F_box, f_box[0])
    result = dict(kind='uniform_profile_box', h=h, R_box=R_box, L_box=L_box, f_box=f_box, F_box=F_box,
                  all_availability_closure_counts=True, moment_log2_integer_upper=moment_bits,
                  global_profile_union_bits=PC_PROFILE_UNION_BITS,
                  no_additional_R_range_union=True, full_code_certificate=False, **checks)
    if not checks['ratios_valid']:
        return dict(**result, certified=False, reason='Required exact corner monotonicity/support check failed')
    R, L, f, F = R_box[0], L_box[1], f_box[0], F_box[0]
    routing = profile_certificate_route(h, R, f, F, cache)
    # Reducing R is valid even through a probability cap.  Reducing F through
    # a capped plateau is not: when corner Q>1 retain the uncapped upper bound.
    routing_bits = routing['probability_capped_bits'] if routing['raw_ratio_at_most_one'] else routing['routebits']
    exponent = _profile_certificate_exponent(moment_bits, R, L, f, F, routing_bits)
    return dict(**result, certified=exponent <= PC_THRESHOLD,
                corner=dict(R=R, L=L, f=f, F=F), routing_cap_used=routing['raw_ratio_at_most_one'],
                exact_uniform_profile_exponent_upper=exponent,
                qualifies_global_profile_class=exponent <= PC_THRESHOLD,
                qualifying_threshold=PC_THRESHOLD, **routing)


def profile_certificate_recompute_moment(iterations=256, mantissa_bits=48, progress=False):
    for name in ('run_certificate_propose', 'run_certificate_verify', 'run_occupancy_matrix'):
        if name not in globals():
            raise RuntimeError('Inject run-occupancy and run-certificate definitions before --certify')
    profile_certificate_constants()
    started = time.monotonic()
    vector, proposal = run_certificate_propose(PC_WIDTH, iterations, mantissa_bits)
    if progress:
        print(dict(stage='profile_moment_witness_proposed', seconds=round(time.monotonic() - started, 3),
                   floating_proposal_only=True), flush=True)
    t = Fraction(PC_Q + 65536, 256 * (PC_Q + 257))
    certificate = run_certificate_verify(vector, PC_T, s=Fraction(75, 64), t=t,
                                         target_lambda_log=5, max_lambda_log=5)
    if not certificate['verified']:
        raise ArithmeticError('Exact row verification did not certify the requested moment multiplier')
    if (certificate['width'] != PC_WIDTH or certificate['updates'] != PC_T
            or certificate['s_numerator'] != 75 or certificate['s_denominator'] != 64
            or Fraction(certificate['t_numerator'], certificate['t_denominator']) != t
            or certificate['occupancy_quartet_numerator'] != 1
            or certificate['occupancy_quartet_denominator'] != 2):
        raise ArithmeticError('Moment certificate parameters do not match the profile theorem')
    return certificate


def profile_certificate_checks():
    constants = profile_certificate_constants()
    logs = recurrences = corners = capped_corners = 0
    for numerator in range(1, 65):
        for denominator in range(1, 65):
            exponent = profile_certificate_ceil_log2_ratio(numerator, denominator)
            assert Fraction(2) ** (exponent - 1) < Fraction(numerator, denominator) <= Fraction(2) ** exponent
            logs += 1
    cache = {}
    for n in range(4, 35):
        for k in range(n + 1):
            assert _profile_certificate_choose(n, k, cache) == comb(n, k)
            recurrences += 1
    for N, selected, base, q in ((24, 8, 4, 16), (24, 8, 4, 256), (16, 8, 4, 256)):
        denominator = comb(N, selected)
        for Flo in range(5, N - 2):
            Fhi = Flo + 1
            for Rlo in range(base, base + N - selected):
                Rhi = Rlo + 1
                flo = 1
                test = profile_certificate_box_ratios(N, selected, base, q, (Rlo, Rhi), (Flo, Fhi), flo)
                if not test['ratios_valid']:
                    continue
                C = Flo - Rlo + base
                corner_Q = Fraction(comb(Flo - flo, C - flo) * comb(N - Flo, selected - C), denominator)
                corner = corner_Q ** 4 * Fraction(2) ** (Flo - flo) / q ** (4 * Rlo)
                for R in range(Rlo, Rhi + 1):
                    for F in range(Flo, Fhi + 1):
                        for f in range(flo, min(test['C_min'], 3) + 1):
                            C = F - R + base
                            Q = Fraction(comb(F - f, C - f) * comb(N - F, selected - C), denominator)
                            actual = Q ** 4 * Fraction(2) ** (F - f) / q ** (4 * R)
                            assert actual <= corner
                            corners += 1
                            # This denominator counts all selected subsets, so Q<=1 here.
                            assert Q <= 1 and corner_Q <= 1
                            capped_corners += 1
    assert (24 * (PC_N + 1) - (117 * (PC_N + 1)) // 512
            > 24 * PC_N - (117 * PC_N) // 512)
    print(dict(integer_ratio_checks=logs, exact_binomial_recurrence_checks=recurrences,
               exact_toy_corner_inequalities=corners, exact_capped_corner_inequalities=capped_corners,
               exact_constants_verified=True, moment_recomputed=False,
               full_code_certificate=False, **constants), flush=True)


def profile_certificate_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--R', type=int, default=900000)
    parser.add_argument('--L', type=int, default=900000)
    parser.add_argument('--f', type=int, default=20000)
    parser.add_argument('--F', type=int, default=1150000)
    parser.add_argument('--half-width', type=int, default=512)
    parser.add_argument('--iterations', type=int, default=256)
    parser.add_argument('--mantissa-bits', type=int, default=48)
    args = parser.parse_args()
    if args.half_width < 0:
        parser.error('--half-width must be nonnegative')
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
            profile_certificate_checks()
            return
        certificate = profile_certificate_recompute_moment(args.iterations, args.mantissa_bits, progress=True)
        moment_bits = certificate['moment_log2_integer_upper']
        print(dict(stage='profile_moment_verified', seconds=round(time.monotonic() - started, 3), **certificate), flush=True)
        cache = {}
        point = profile_certificate_point(args.h, args.R, args.L, args.f, args.F, moment_bits, cache)
        print(dict(seconds=round(time.monotonic() - started, 3), **point), flush=True)
        radius = args.half_width
        box = profile_certificate_box(args.h, (args.R - radius, args.R + radius),
                                       (args.L - radius, args.L + radius),
                                       (args.f - radius, args.f + radius),
                                       (args.F - radius, args.F + radius), moment_bits, cache)
        print(dict(seconds=round(time.monotonic() - started, 3), **box), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    profile_certificate_main()
