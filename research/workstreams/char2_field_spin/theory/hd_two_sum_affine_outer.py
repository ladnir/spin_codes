"""Exact affine two-sum algebra and one bounded sparse-cohort FLOAT screen.

H=I+11^T+uu^T with u=(8,...,15) in AES GF256 is a NEW outer ensemble.
Every singleton has full output support. Every fixed input pair has the
same law: full support with probability247/255 and each one-hole support
with probability1/255. Independent output diagonals restore independent
nonzero amplitudes conditional on support. Fixed row sets are not resampled.

The positive ordered hole recurrence removes the old pair-multiplicity
resource relaxation. It still uses established source-conditioned regional
upper matrices and one common potential. No old certificate transfers.
The optional CLI runs exactly two frozen witnesses at r17/source(8,0), all
18 sparse cohort splits. No fit, native code, benchmark or certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from gf8_packet_maps import MUL, rank, xor_sum
from hd_two_sum_outer import silent_moment_count


LABELS = tuple(range(8, 16))
CONSTRUCTION = 'two-sum-affine-gf256-h4-alternating-swapped-pairs-k20-rate-half'
WITNESSES = (
    ('inherited', F('996248217/1000000000'), F('-127/2')),
    ('sparse-mixed-refit', F('4487948476156805/4503599627370496'), F('-8814329/131072')),
)


@lru_cache(None)
def matrix():
    return tuple(tuple(int(i == j) ^ 1 ^ MUL[LABELS[i]][LABELS[j]] for j in range(8))
                 for i in range(8))


def _inputs(values):
    if len(values) != 8 or any(type(value) is not int or not 0 <= value < 256 for value in values):
        raise ValueError('Eight AES GF256 bytes required')


def direct(values):
    _inputs(values)
    s = xor_sum(values)
    t = xor_sum(MUL[label][value] for label, value in zip(LABELS, values))
    return tuple(value ^ s ^ MUL[label][t] for label, value in zip(LABELS, values))


def circuit(values):
    """Six fixed products and30 binary XORs; the native circuit is not edited."""
    _inputs(values)
    x0, x1, x2, x3, x4, x5, x6, x7 = values
    p01, p23, p45, p67 = x0 ^ x1, x2 ^ x3, x4 ^ x5, x6 ^ x7
    b = p23 ^ p67
    s = (p01 ^ p45) ^ b
    c = p45 ^ p67
    a = (x1 ^ x3) ^ (x5 ^ x7)
    t = a ^ MUL[2][b] ^ MUL[4][c] ^ MUL[8][s]
    sprime = s ^ MUL[8][t]
    t2, t4 = MUL[2][t], MUL[4][t]
    c0, c1, c2 = sprime, sprime ^ t, sprime ^ t2
    c3, c4, c5 = c1 ^ t2, sprime ^ t4, c1 ^ t4
    c6, c7 = c2 ^ t4, c3 ^ t4
    return tuple(x ^ correction for x, correction in zip(values, (c0, c1, c2, c3, c4, c5, c6, c7)))


@lru_cache(None)
def uniform_mask_counts(active):
    """Exact uniform-coordinate input counts, zeros allowed; lazy full kernel."""
    if type(active) is not int or not 0 <= active < 256:
        raise ValueError('Eight-bit active coordinate mask required')
    columns = [i for i in range(8) if active >> i & 1]
    h = matrix()
    counts = [256**(len(columns)-rank([[h[i][j] for j in columns]
              for i in range(8) if not target >> i & 1])) for target in range(256)]
    for bit in range(8):
        for target in range(256):
            if target >> bit & 1:
                counts[target] -= counts[target ^ (1 << bit)]
    if min(counts) < 0 or sum(counts) != 256**len(columns):
        raise ArithmeticError('Exact output-support inversion failed')
    return tuple(counts)


@lru_cache(None)
def nonzero_mask_counts(active):
    """Exact kernel with every selected input nonzero, via input inclusion-exclusion."""
    if type(active) is not int or not 0 <= active < 256:
        raise ValueError('Eight-bit active coordinate mask required')
    result = [0]*256
    subset = active
    while True:
        sign = -1 if (active.bit_count()-subset.bit_count()) % 2 else 1
        for target, count in enumerate(uniform_mask_counts(subset)):
            result[target] += sign*count
        if not subset:
            break
        subset = (subset-1) & active
    if min(result) < 0 or sum(result) != 255**active.bit_count():
        raise ArithmeticError('Exact input-support inversion failed')
    return tuple(result)


def sparse_mask_law(arity):
    """Exact conditional support probabilities for any fixed support of size<=2."""
    if type(arity) is not int or not 0 <= arity <= 2:
        raise ValueError('This uniform sparse law covers only arity0,1,2')
    if arity == 0:
        return {0: F(1)}
    if arity == 1:
        return {255: F(1)}
    return {255: F(247, 255), **{255 ^ (1 << i): F(1, 255) for i in range(8)}}


def column_action_exact(regions, potential, single, double, q=256):
    """Positive binomial EGF recurrence for ordered regional matrices.

    For t=single+double and k labelled holes distributed across suffix rows,
    V_i[k]=sum_(h<=k) C(k,h) R_i(t-h)V_(i+1)[k-h]. Start V_L[0]=v,
    V_L[k>0]=0. Mix sum_k C(d,k)*(q-1-L)^(d-k)V_0[k]/(q-1)^d.
    This is exact for the supplied matrices; no regional resource maximum
    or potential quotient is taken here. Noncommuting order is preserved.
    """
    length, t = len(regions), single+double
    if (type(single) is not int or type(double) is not int or min(single, double) < 0
            or not 1 <= length < q or len(potential) != 2
            or any(len(table) <= t for table in regions)):
        raise ValueError('Complete ordered tables and valid sparse category required')
    tables = tuple(tuple(tuple(F(entry) for entry in row) for row in table) for table in regions)
    vector = tuple(F(value) for value in potential)
    if any(value <= 0 for value in vector) or any(len(row) != 4 or min(row) < 0 for table in tables for row in table):
        raise ValueError('Positive potential and nonnegative flat-four matrices required')
    values = [vector]+[(F(0), F(0))]*double
    for table in reversed(tables):
        following = []
        for k in range(double+1):
            first, second = F(0), F(0)
            for h in range(k+1):
                a, b, c, d = table[t-h]
                left, right = values[k-h]
                factor = comb(k, h)
                first += factor*(a*left+b*right)
                second += factor*(c*left+d*right)
            following.append((first, second))
        values = following
    return tuple(sum((F(comb(double, k)*(q-1-length)**(double-k), (q-1)**double)*values[k][j]
                      for k in range(double+1)), F(0)) for j in range(2))


def column_action_float(regions, log_potential, single, double, q=256):
    """Log2 implementation of the same positive ordered recurrence; FLOAT only."""
    tables = np.asarray(regions, dtype=float)
    t = single+double
    if (tables.ndim != 4 or tables.shape[2:] != (2, 2) or not 1 <= len(tables) < q-1
            or type(single) is not int or type(double) is not int or min(single, double) < 0
            or t >= tables.shape[1] or not np.isfinite(log_potential)
            or np.any(np.isnan(tables)) or np.any(np.isposinf(tables))):
        raise ValueError('Complete finite-or-zero log-matrix tables and valid category required')
    log_binomials = [np.array([log2(comb(k, h)) for k in range(h, double+1)]) for h in range(double+1)]
    values = np.full((double+1, 2), -np.inf)
    values[0] = [0., log_potential]
    for table in reversed(tables):
        following = np.full_like(values, -np.inf)
        for h in range(double+1):
            matrix, previous = table[t-h], values[:double+1-h]
            action = np.logaddexp2(matrix[None, :, 0]+previous[:, 0, None],
                                  matrix[None, :, 1]+previous[:, 1, None])
            np.logaddexp2(following[h:], action+log_binomials[h][:, None], out=following[h:])
        values = following
    probabilities = np.array([log2(comb(double, k))+(double-k)*log2(q-1-len(tables))
                              -double*log2(q-1) for k in range(double+1)])
    return np.logaddexp2.reduce(values+probabilities[:, None], axis=0)


def sparse_column_coefficients(caps, n1, n2, q=256):
    """Exact positive coefficients, including(q-1)^X; row-set choices are outside.

    caps[(single,double)] bounds the complete chronological column action
    under one common potential. Raise this polynomial to32, retain degree
    >=17*(n1+2*n2), then divide by(q-1)^(16*(n1+2*n2)).
    """
    if any(type(n) is not int or n < 0 for n in (n1, n2)) or q < 2:
        raise ValueError('Nonnegative sparse cohort sizes required')
    result = [F(0)]*(n1+2*n2+1)
    for a in range(n1+1):
        for s in range(n2+1):
            for d in range(n2-s+1):
                cap = F(caps[a+s, d])
                if cap < 0:
                    raise ValueError('Nonnegative complete category caps required')
                degree = a+s+2*d
                result[degree] += comb(n1, a)*comb(n2, s)*comb(n2-s, d)*2**s*(q-1)**degree*cap
    return tuple(result)


def screen(progress=None):
    """Only r17, source(8,0), all18 sparse cohorts and two existing literals."""
    from hd_xor_group_geometry import matching_joint_counts
    from hd_xor_joint_columns import categories
    from hd_xor_mixed_gate import ordered_region_tables, sparse_log_count
    from hd_xor_multigroup import LIMIT
    from hd_xor_pair_route import pair_regions_screen

    started = perf_counter()
    source_count = matching_joint_counts(8192, 17)[8, 0]
    banks, selected = [], {}
    for label, z, logp in WITNESSES:
        bank_started = perf_counter()
        raw = pair_regions_screen(float(z), 17)
        regions = ordered_region_tables(raw, 17, 8, 0)
        caps = {}
        for single, double in categories(17):
            vector = column_action_float(regions, float(logp), single, double)
            caps[single, double] = float(np.max(vector-[0., float(logp)]))
        entries = []
        for n1 in range(18):
            n2 = 17-n1
            value = (log2(source_count)+log2(comb(17, n1))+3*n1+n2*log2(28)
                     -log2(255)-LIMIT*log2(float(z))-float(logp)+sparse_log_count(caps, n1, n2))
            row = dict(cohort_counts=[n1, n2, 0], witness_label=label, mark=str(z), log2_potential=str(logp),
                source_supports=source_count, source_colorings=comb(17, n1),
                singleton_row_choices=8**n1, pair_row_choices=28**n2,
                sparse_support_floor=17*(n1+2*n2), sparse_normalizer_exponent=16*(n1+2*n2),
                failure_log2_screen=value)
            entries.append(row)
            if n1 not in selected or value < selected[n1]['failure_log2_screen']:
                selected[n1] = row
            if progress:
                progress(dict(stage='cohort', **row))
        banks.append(dict(witness_label=label, mark=str(z), log2_potential=str(logp), entries=entries,
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
            elapsed_seconds=perf_counter()-bank_started))
    chosen = [selected[n1] for n1 in range(18)]
    return dict(construction=CONSTRUCTION, field_polynomial='0x11b', labels=list(LABELS),
        field_order=256, forbidden_weight=LIMIT, active_groups=17, source_matching_counts=[8, 0],
        source_classes=1, all_source_classes=81, total_columns=32, complete_source_coverage=False,
        complete_sparse_cohorts_for_selected_source=True, cohort_count=18,
        pair_swaps=True, fixed_rows_across_columns=True, fixed_source_geometry=True,
        local_mask_mixture_exact=True, regional_matrices_are_upper_bounds=True,
        potential_quotient_after_hole_mixture=True, ordered_regional_matrices=True,
        pair_multiplicity_strata_used=False, resource_maximum_used=False,
        witness_banks=2, maximum_witness_banks=2, fitted_new_witnesses=0,
        projective_normalizations=1, endpoint_potential_factors=1,
        banks=banks, selected=chosen,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in chosen])),
        failed_cohorts_at_40_bits=[row['cohort_counts'] for row in chosen if row['failure_log2_screen'] >= -40],
        directed=False, old_certificates_transferred=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', action='store_true', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded research artifact')
    result = screen(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('banks', 'selected')}, allow_nan=False))
