"""Exact death slices of the fixed hard ACTUAL-state profile.

Profile: h=4096,R=Ls=786432,f=4,F=1048576. The existing actual-state
census complements remove none of its physical D interval [0,262144].
The new comparisons prove, for this fixed profile:
  D<=32768: probability <2^-14587;
  D>=196608: probability <2^-26713.
The optional larger low slice D<=34569 has exponent -179.

Both counts already sum their full death interval: there is NO D union.
These are slices of one profile, not a whole-profile or full-code proof.
The conventional 96-bit global-profile union is reported separately;
numeric bounds for this one profile do not assert coverage elsewhere.

--checks is a small self-contained integer/rational derivation check. Its
LOW slice invokes the previously proved 22323-bit run/occupancy theorem;
it does not pretend to reverify that matrix witness. Reference SHA256s
identify the frozen source supporting that theorem, not a saved witness.
The HIGH slice and routing comparison do not require the moment theorem.

--certify additionally proposes and verifies a fresh positive integer
vector. Inject quartet_run_occupancy.py and quartet_run_certificate.py
definitions with main blocks removed before using that mode. All 513
rows are rechecked with integers; the floating vector is only a proposal.
No witness files, dumps, local imports or subprocesses. Run jobs serially.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
from math import comb
import os
import shutil
import threading
import time


HDS_N = 1 << 21
HDS_Q = 1 << 32
HDS_WIDTH = 512
HDS_UPDATES = 4096
HDS_W = 943718
HDS_H = 4096
HDS_R = HDS_L = 786432
HDS_F = 1048576
HDS_BIRTHS = 4
HDS_MOMENT_BITS = 22323
HDS_UNION_BITS = 96
HDS_REFERENCE_SHA256 = {
    'quartet_run_occupancy.py': '3e066747f42a56fa25bda16b42c8ff53677419d3ecaf550a9985b5d46d6c285e',
    'quartet_run_certificate.py': 'ae08f37541905bc68eaf63e824f116ed7c546f686a3c9cb747044c0a706c5dff',
}


def _hds_census_constant(j):
    if j == 0:
        return HDS_N
    if j == 3:
        return 361984
    denominator = 1 << (j + 1)
    return HDS_WIDTH * j + (3 * HDS_N + denominator - 1) // denominator


def hard_death_slices_checks():
    """Check every finite comparison used by the stated analytic slices."""
    N, q, L, R, F, f, W = HDS_N, HDS_Q, HDS_L, HDS_R, HDS_F, HDS_BIRTHS, HDS_W
    assert N == HDS_WIDTH * HDS_UPDATES and R == L
    assert F == N // 2 and L == 3 * N // 8
    G = F - f
    J, k = (N - HDS_WIDTH) // 4, G // 4
    assert G % 4 == 0 and (J, k) == (524160, 262143)
    assert 0 <= k <= J and k <= L <= G + HDS_WIDTH
    C_H, C_out = F - R + 128 * HDS_H, R - 128 * HDS_H
    assert C_H == 3 * F // 4 and C_out == F // 4
    assert 129 * HDS_H <= C_H <= 256 * HDS_H
    assert 0 <= C_out <= 128 * (8192 - HDS_H)
    assert f <= C_H and 256 * HDS_H - C_H == F // 4
    d_min, d_max = max(0, f - HDS_WIDTH), min(F - L, W // 2)
    assert (d_min, d_max) == (0, 262144)

    # Every implemented rank/closure exclusion fails throughout the interval.
    census_checks = 0
    assert R - L < 65537
    for j in range(17):
        c = _hds_census_constant(j)
        rank_max = 32 * (R - L) - 2 * j * d_min
        assert rank_max < c + 45
        coefficient = 16 - 2 * j
        closure_max = 32 * R - 16 * L + coefficient * (d_max if coefficient >= 0 else d_min)
        assert closure_max < 16 * W + c + 92
        census_checks += 2

    # Q<=C(F,F/4)^2/C(8192,4096)^255. The forced-birth restriction only
    # reduces the first numerator. Binomial evaluation at 1/3 gives its bound.
    assert F % 8 == 0 and 3 ** 12 > 1 << 19
    assert 8193 < 1 << 14
    central = comb(8192, 4096)
    assert 8193 * central > 1 << 8192
    assert central > 1 << 8178
    route_numerator_bits = 13 * F // 8
    route_denominator_bits = 8178 * 255
    route_bits = route_numerator_bits - route_denominator_bits
    assert (route_numerator_bits, route_denominator_bits, route_bits) == (1703936, 2085390, -381454)

    # The rank and exact optional-coefficient bounds both use this correction.
    correction = Fraction(1024, 989) / (1 - Fraction(N, q))
    assert correction == Fraction(1024, 989) * Fraction(2048, 2047) < 2
    rank_field_bits = -31 + 32 * (L - R)
    assert rank_field_bits == -31

    # The LOW slice relies on the published exact 513-row moment theorem.
    s = Fraction(75, 64)
    t = Fraction(q + 65536, 256 * (q + 257))
    assert Fraction(1, 256) <= t < 1
    assert 75 ** 1024 > 1 << (6 * 1024 + 234)
    assert s ** 1024 > 1 << 234
    assert L % 1024 == 0 and L % 512 == 0
    s_discount_bits = 117 * L // 512
    assert s_discount_bits == 179712
    assert HDS_MOMENT_BITS == 5 * HDS_UPDATES + 1843
    low_count_constant = HDS_MOMENT_BITS - s_discount_bits + k
    low_probability_constant = low_count_constant + route_bits + rank_field_bits
    assert (low_count_constant, low_probability_constant) == (104754, -276731)
    low_cut = 32768
    low_bits = low_probability_constant + 8 * low_cut
    larger_low_cut = (-172 - low_probability_constant) // 8
    larger_low_bits = low_probability_constant + 8 * larger_low_cut
    assert (low_bits, larger_low_cut, larger_low_bits) == (-14587, 34569, -179)
    assert low_probability_constant + 8 * (larger_low_cut + 1) > -172

    # HIGH geometry: coefficient bound at s=3, with the free terminal 4^512.
    assert 255 < 1 << 8 and L % 12 == 0
    assert 3 ** 12 > 1 << 19
    geometry_bits = J + 8 * k + 2 * HDS_WIDTH - 19 * (L // 12)
    assert geometry_bits == 1377144
    death_threshold = W - L
    high_cut = 196608
    delta = high_cut - death_threshold
    n_optional = L - high_cut
    assert (death_threshold, delta, n_optional) == (157286, 39322, 589824)
    assert death_threshold < high_cut <= d_max
    assert n_optional <= 20 * delta
    assert 17 ** 10 < 2 * 16 ** 10
    # At every D>=high_cut, n=L-D decreases and Delta increases. Therefore
    # C(n,Delta)<= (17/16)^n*16^Delta < 2^(6Delta), from the binomial theorem.
    high_field_bits = -31 - (32 - 6) * delta
    high_bits = geometry_bits + route_bits + high_field_bits
    assert (high_field_bits, high_bits) == (-1022403, -26713)

    # An optional endpoint-series interface at D=200000: no death union.
    endpoint = 200000
    run_index = endpoint + HDS_WIDTH
    zeros = N - L + HDS_WIDTH
    endpoint_delta = endpoint - death_threshold
    endpoint_n = L - endpoint
    assert L - run_index < 3 * run_index
    assert zeros - run_index < 6 * run_index
    assert endpoint_n < 14 * (endpoint_delta + 1)
    assert 18 * 14 < 1 << 8
    # The adjacent combined run/field ratio is <2^8/q=2^-24.

    assert low_bits < -172 and high_bits < -172 and larger_low_bits < -172
    assert 8192 * (1 << 20) * N ** 3 == 1 << HDS_UNION_BITS
    return dict(statistics='actual_state_support', profile=(HDS_H, R, L, f, F),
                physical_death_interval=(d_min, d_max), census_exclusions_checked=census_checks,
                census_complement_death_interval=(d_min, d_max),
                routing_exponent_upper=route_bits,
                low_slice=(d_min, low_cut), low_slice_exponent_upper=low_bits,
                optional_larger_low_slice=(d_min, larger_low_cut),
                optional_larger_low_exponent_upper=larger_low_bits,
                high_slice=(high_cut, d_max), high_slice_exponent_upper=high_bits,
                middle_interval_still_open=(larger_low_cut + 1, high_cut - 1),
                global_profile_union_bits=HDS_UNION_BITS,
                combined_tails_after_optional_profile_allowance=max(larger_low_bits, high_bits) + 1 + HDS_UNION_BITS,
                no_additional_death_union=True,
                prior_moment_log2_bound=HDS_MOMENT_BITS, prior_moment_reverified=False,
                low_slice_conditional_on_prior_moment_theorem=True,
                high_slice_independent_of_moment_theorem=True,
                reference_source_sha256=HDS_REFERENCE_SHA256,
                source_hash_is_not_a_witness_certificate=True,
                whole_profile_certificate=False, full_code_certificate=False)


def hard_death_slices_reverify():
    """Optional fresh exact support for the invoked 22323-bit moment theorem."""
    for name in ('run_occupancy_matrix', 'run_certificate_propose', 'run_certificate_verify'):
        if name not in globals():
            raise RuntimeError('Inject run/occupancy and run-certificate definitions; missing ' + name)
    vector, proposal = run_certificate_propose(512, iterations=256, mantissa_bits=48)
    q = HDS_Q
    s = Fraction(75, 64)
    t = Fraction(q + 65536, 256 * (q + 257))
    certificate = run_certificate_verify(vector, updates=4096, s=s, t=t,
                                         target_lambda_log=5, max_lambda_log=5)
    if not (certificate['verified'] and certificate['exact_rows'] == 513
            and certificate['width'] == 512 and certificate['updates'] == 4096
            and certificate['lambda_log2'] == 5
            and certificate['moment_log2_integer_upper'] <= HDS_MOMENT_BITS
            and certificate['s_numerator'] == 75 and certificate['s_denominator'] == 64
            and Fraction(certificate['t_numerator'], certificate['t_denominator']) == t
            and certificate['occupancy_quartet_numerator'] == 1
            and certificate['occupancy_quartet_denominator'] == 2):
        raise ArithmeticError('The fresh integer witness does not establish the invoked moment theorem')
    digest = sha256(','.join(str(value) for value in vector).encode('ascii')).hexdigest()
    return dict(moment_certificate=certificate, witness_sha256=digest,
                witness_hash_encoding='comma-separated positive decimal integers, no trailing comma',
                floating_proposal_only=True, prior_moment_reverified=True)


def hard_death_slices_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        started = time.monotonic()
        result = hard_death_slices_checks()
        if args.certify:
            result.update(hard_death_slices_reverify())
        print(dict(**result, seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    hard_death_slices_main()
