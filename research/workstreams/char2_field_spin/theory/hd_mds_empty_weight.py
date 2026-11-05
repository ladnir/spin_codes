"""Exact empty-packet weight law and a bounded FLOAT power diagnostic.

H is the INNER16-point interpolation map, not the rank-five outer H*.
For fresh independent nonzero diagonals D1,D2, its superregularity makes
D2 H D1 x uniform within each output-weight class, for every fixed x of
a given weight. Weighted histories therefore close on17 weight classes.

For X=0 the emission precedes refresh, and K[a,b]=e_a(z) P[a,b].
e_a is computed from the actual four RS[16,4] coefficient-basis maps;
no selected coefficient-column MDS assumption is used. The zero state
is absorbing, and nonzero states cannot return to zero in empty packets.

This diagnostic builds no nonempty-packet or regional support provider.
The integer counts and rational kernels are exact; eigenvalues/powers in
screen() are FLOAT proposals, not directed probability certificates.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_grouped_quiet import quiet_cap, weight_transition
from hd_xor_state_local import batch_ranks, checked_emission, current_emission


def local_positive_counts(emission, q=256, *, check_time=lambda: None):
    """Counts N[j,d]=sum_|Z|d #{wt(Q)=j : (A Q)_Z=0}, j>0.

    The full map is MDS, so nonzero Q gives at most h-1 output zeros.
    For every selected Z, input-coordinate inclusion-exclusion counts
    exact input weight using actual ranks of A[Z,J]. The empty state is
    recorded separately, because it has n rather than h-1 output zeros.
    """
    emission = tuple(map(tuple, emission))
    n, h = checked_emission(emission, q)
    result = [[0]*h for _ in range(h+1)]
    for zeros in range(h):
        check_time()
        choices = tuple(combinations(range(n), zeros))
        matrices = np.asarray([[emission[i] for i in positions] for positions in choices],
                              dtype=np.uint8).reshape(len(choices), zeros, h)
        for mask in range(1 << h):
            columns = tuple(i for i in range(h) if mask >> i & 1)
            dimension = len(columns)
            ranks = (batch_ranks(matrices[:, :, columns], q)
                     if zeros and dimension else np.zeros(len(choices), dtype=np.int64))
            solutions = sum(q**(dimension-int(rank)) for rank in ranks)
            for weight in range(max(1, dimension), h+1):
                result[weight][zeros] += (-1)**(weight-dimension)*comb(h-dimension, weight-dimension)*solutions
    if any(value < 0 for row in result for value in row):
        raise ArithmeticError('Negative exact positive zero-set count')
    if any(result[j][0] != comb(h, j)*(q-1)**j for j in range(1, h+1)):
        raise ArithmeticError('Input-weight marginal does not normalize')
    return tuple(map(tuple, result))


def local_marked_counts(emission, z, q=256, *, check_time=lambda: None):
    z = F(z)
    if not 0 < z <= 1:
        raise ValueError('Exact mark in(0,1] required')
    n, h = len(emission), len(emission[0])
    counts = local_positive_counts(emission, q, check_time=check_time)
    return (F(1),)+tuple(sum((F(counts[j][d])*z**(n-d)*(1-z)**d
                             for d in range(h)), F(0)) for j in range(1, h+1))


def emission_by_weight(z, *, q=256, blocks=4, emission=None, check_time=lambda: None):
    """Exact marked mean in each uniformly supported nonzero-amplitude class."""
    emission = current_emission() if emission is None else tuple(map(tuple, emission))
    if type(blocks) is not int or not 1 <= blocks <= 4:
        raise ValueError('One to four equal blocks required')
    h = len(emission[0])
    local = local_marked_counts(emission, z, q, check_time=check_time)
    coefficients = [F(1)]
    for _ in range(blocks):
        check_time()
        following = [F(0)]*(len(coefficients)+h)
        for i, left in enumerate(coefficients):
            for j, right in enumerate(local):
                following[i+j] += left*right
        coefficients = following
    state = blocks*h
    return tuple(value/(comb(state, j)*(q-1)**j) for j, value in enumerate(coefficients))


def empty_kernel(z, *, q=256, blocks=4, emission=None, check_time=lambda: None):
    emission = current_emission() if emission is None else tuple(map(tuple, emission))
    state = blocks*len(emission[0])
    means = emission_by_weight(z, q=q, blocks=blocks, emission=emission, check_time=check_time)
    transition = weight_transition(q, state)
    zero = F(0)
    kernel = ((F(1),)+(zero,)*state,)+tuple((zero,)+tuple(means[a+1]*value for value in row)
                                           for a, row in enumerate(transition))
    return kernel, means, transition


def flat_cap_from_kernel(means, transition, q):
    """Exact one-packet cap in the UNCHANGED uniform-live density cone.

    Inverse refresh is also two NZ diagonals around a superregular map.
    At target weight b, its marked expectation is sum_a P[b,a] e_a.
    Detailed balance gives the same quantity as (pi K)[b]/pi[b].
    """
    state = len(transition)
    if (len(means) != state+1 or means[0] != 1
            or any(len(row) != state for row in transition)):
        raise ValueError('Matching exact weight means and live transition required')
    pi = tuple(F(comb(state, a)*(q-1)**a, q**state-1) for a in range(1, state+1))
    inverse_values = tuple(sum((p*value for p, value in zip(row, means[1:])), F(0))
                           for row in transition)
    for b, value in enumerate(inverse_values):
        forward = sum((pi[a]*means[a+1]*transition[a][b] for a in range(state)), F(0))/pi[b]
        if value != forward:
            raise ArithmeticError('Forward density and inverse-refresh expectations disagree')
    return max(inverse_values), inverse_values, pi


def exact_scalar_receipt(z=F(9927, 10000), *, max_seconds=10):
    """Exact rational empty11 plus directed logs; no regional provider.

    The exact rational inequality is valid pointwise at every live output
    state. All nonempty packet entries can stay unchanged. The128th power
    below covers only a completely empty region, not other histograms.
    """
    from flint import arb, ctx
    if not 0 < max_seconds <= 10:
        raise ValueError('Explicit budget at most10 seconds required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Exact empty scalar budget exhausted')
    kernel, means, transition = empty_kernel(z, check_time=check_time)
    scalar, targets, pi = flat_cap_from_kernel(means, transition, 256)
    beta = F(256**16-1, 255**16)
    old = min(beta*sum((p*e for p, e in zip(pi, means[1:])), F(0)), quiet_cap(z))
    if not 0 < scalar <= old <= 1:
        raise ArithmeticError('Exact empty scalar does not improve the old cap')
    def rational(value):
        return arb(value.numerator)/value.denominator
    def endpoint(value, upper):
        value = value if value.is_exact() else (value.upper() if upper else value.lower())
        exact = value.fmpq()
        return F(int(exact.p), int(exact.q))
    directed = []
    denominator = 10**12
    for precision in (256, 384):
        check_time()
        with ctx.workprec(precision):
            new_log = 128*rational(scalar).log()/arb(2).log()
            old_log = 128*rational(old).log()/arb(2).log()
            upper = endpoint(new_log, True)
            ceiling = F(-((-upper.numerator*denominator)//upper.denominator), denominator)
            gain = endpoint(old_log-new_log, False)
            floor = F((gain.numerator*denominator)//gain.denominator, denominator)
            directed.append(dict(precision_bits=precision, new_empty_region_log2_ceiling=str(ceiling),
                                 improvement_bits_floor=str(floor)))
    check_time()
    return dict(mark=str(z), field_order=256, state=16, packet_width=64,
        emission='four actual normalized RS[16,4] blocks; coefficient columns not assumed MDS',
        refresh='fresh independent NZ D2 H16 D1; superregular interpolation H16',
        exact_empty11=str(scalar), old_empty11=str(old),
        maximizing_target_weights=[b+1 for b, value in enumerate(targets) if value == scalar],
        all_16_forward_inverse_ratios_equal=True, unchanged_live_reference='uniform on GF256^16 excluding0',
        nonempty_packet_caps_unchanged=True, directed_comparisons=directed,
        scope='empty-packet11 replacement;128-packet comparison only for a completely empty region',
        regional_provider_built=False, whole_code_certificate=False,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-started)


def screen(z=F(9927, 10000), *, max_seconds=30):
    if not 0 < max_seconds <= 30:
        raise ValueError('A budget of at most30 seconds is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Empty-weight diagnostic budget exhausted')
    kernel, means, transition = empty_kernel(z, check_time=check_time)
    state, q, packets = 16, 256, 128
    pi = tuple(F(comb(state, a)*(q-1)**a, q**state-1) for a in range(1, state+1))
    # Independently check the uniform-nonzero stationary law exactly.
    for b in range(state):
        if sum(pi[a]*transition[a][b] for a in range(state)) != pi[b]:
            raise ArithmeticError('Exact stationary law failed')
    live = np.array([[float(value) for value in row[1:]] for row in kernel[1:]])
    distribution = np.array(list(map(float, pi)))
    matrix = np.linalg.matrix_power(live, packets)
    finals = distribution @ matrix
    total = float(finals.sum())
    # This collapses the17-component result back into the SAME normalized
    # uniform-nonzero cone. It compares densities, not just total mass.
    density = float(np.max(finals/distribution))
    beta = F(q**state-1, (q-1)**state)
    average = sum(weight*mean for weight, mean in zip(pi, means[1:]))
    old = min(beta*average, quiet_cap(z))
    old_region = packets*log2(old)
    spectral = float(np.max(np.abs(np.linalg.eigvals(live))))
    check_time()
    return dict(mark=str(z), scope='128 consecutive empty packets; original one-state inner only',
        refresh='D2 H16 D1; superregular interpolation H16; fresh independent NZ diagonals',
        state_order=list(range(17)), packet_marked_means=list(map(float, means)),
        exact_transition_stationarity_checked=True, original_live_reference='uniform on GF256^16 excluding0',
        old_packet_log2=log2(old), old_regional_R11_log2=old_region,
        exact_weight_kernel_total_mass_log2=log2(total),
        exact_weight_kernel_uniform_density_log2=log2(density),
        density_gain_bits=old_region-log2(density), worst_endpoint_weight=1+int(np.argmax(finals/distribution)),
        worst_initial_weight_mass_log2=log2(float(np.max(matrix.sum(axis=1)))),
        live_spectral_radius_log2=log2(spectral), live_spectral_128_log2=packets*log2(spectral),
        packet_beta_log2=log2(beta), repeated_beta_128_log2=packets*log2(beta),
        numerical_method='FLOAT positive16x16 powers/eigenvalues of exact rational kernel',
        nonempty_packet_kernel_built=False, regional_support_provider_built=False,
        directed=False, whole_code_certificate=False, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exact-scalar', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        parser.error('Do not overwrite a frozen empty-weight result')
    answer = exact_scalar_receipt() if args.exact_scalar else screen()
    if args.output is not None:
        with args.output.open('x', encoding='utf-8') as handle:
            json.dump(answer, handle, indent=2, allow_nan=False)
            handle.write('\n')
    print(json.dumps(answer, indent=2, allow_nan=False))
