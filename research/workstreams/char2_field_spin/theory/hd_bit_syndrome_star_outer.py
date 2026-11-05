"""Exact local algebra for three star-quadratic additions to H*.

This helper is not a regional provider or a distance certificate.  The three
new syndromes are a0*a1, a0*a2, a0*a3 on the sixteen binary row labels.
The old construction, its samplers, and all old receipts remain unchanged.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb
import json

import hd_bit_syndrome_outer as original
from gf8_packet_maps import MUL, determinant, matmul, matvec

CONSTRUCTION = 'G16-RS16x8-D2-star-quadratic-syndrome-D1-research'
BASIS = original.BASIS
OFFSET = original.OFFSET
LABELS = original.LABELS
QUADRATIC_COEFFICIENTS = {
    'identity': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    'coupled': ((0, 1, 0), (1, 0, 1), (0, 1, 1)),
}


def quadratic_coefficients(variant='identity'):
    if variant not in QUADRATIC_COEFFICIENTS:
        raise ValueError('variant must be identity or coupled')
    return QUADRATIC_COEFFICIENTS[variant]


def construction_identity(variant='identity'):
    quadratic_coefficients(variant)
    return CONSTRUCTION + '-' + variant


def _mask(mask):
    if type(mask) is not int or not 0 <= mask < (1 << 16):
        raise ValueError('expected a sixteen-label support mask')
    return mask


def syndrome_rows():
    """Rows of [1,a0,a1,a2,a3,a0*a1,a0*a2,a0*a3]."""
    return tuple((1,) + tuple(i >> j & 1 for j in range(4))
                 + tuple((i & 1) * (i >> j & 1) for j in (1, 2, 3))
                 for i in range(16))


def coefficient_matrix(variant='identity'):
    b = original.coefficient_matrix()
    c = quadratic_coefficients(variant)
    return tuple(tuple(b[i][j] if i < 5 and j < 5
                       else c[i - 5][j - 5] if i >= 5 and j >= 5 else 0
                       for j in range(8)) for i in range(8))


@lru_cache(maxsize=2)
def matrix(variant='identity'):
    h = original.matrix()
    v = syndrome_rows()
    c = quadratic_coefficients(variant)
    return tuple(tuple(h[i][j] ^ (sum(v[i][5 + k] * c[k][ell] * v[j][5 + ell]
                                      for k in range(3) for ell in range(3)) & 1)
                       for j in range(16)) for i in range(16))


def _values(values):
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < 256
                                for x in values):
        raise ValueError('expected sixteen AES-field bytes')
    return values


def direct(values, variant='identity'):
    return matvec(matrix(variant), _values(values))


def _circuit(x, add, mul, variant='identity'):
    quadratic_coefficients(variant)
    # Preserve the old pair/quarter sums.  The odd-input pair sums are also
    # needed for the old a0 syndrome; naming them introduces no new XOR.
    p0, p1, p2, p3 = add(x[0], x[1]), add(x[2], x[3]), add(x[4], x[5]), add(x[6], x[7])
    p4, p5, p6, p7 = add(x[8], x[9]), add(x[10], x[11]), add(x[12], x[13]), add(x[14], x[15])
    o0, o1, o2, o3 = add(x[1], x[3]), add(x[5], x[7]), add(x[9], x[11]), add(x[13], x[15])
    a0 = add(add(o0, o1), add(o2, o3))
    a1 = add(add(p1, p3), add(p5, p7))
    q0, q1, q2, q3 = add(p0, p1), add(p2, p3), add(p4, p5), add(p6, p7)
    a2, a3 = add(q1, q3), add(q2, q3)
    s = add(add(q0, q1), a3)
    t = add(a0, add(add(mul(s, OFFSET), mul(a1, BASIS[1])),
                       add(mul(a2, BASIS[2]), mul(a3, BASIS[3]))))
    common = add(s, mul(t, OFFSET))
    delta0 = add(t, a0)
    d1, d2, d3 = add(mul(t, BASIS[1]), a1), add(mul(t, BASIS[2]), a2), add(mul(t, BASIS[3]), a3)

    # Five new syndrome XORs, then three new output-delta XORs.
    b1 = add(add(x[3], x[7]), add(x[11], x[15]))
    b2, b3 = add(o1, o3), add(o2, o3)
    if variant == 'coupled':
        # C*b=(b2,b1+b3,b2+b3): two more XORs, no field product.
        b1, b2, b3 = b2, add(b1, b3), add(b2, b3)
    odd_delta = (add(d1, b1), add(d2, b2), add(d3, b3))
    even_delta = (d1, d2, d3)

    # Two affine eight-point evaluations use 7 + (1+7) XORs, equal to the
    # old sixteen-point correction's 15 XORs.  Label order is unchanged.
    even = [common] * 8
    odd = [add(common, delta0)] * 8
    for i in range(1, 8):
        bit = i & -i
        index = bit.bit_length() - 1
        even[i] = add(even[i ^ bit], even_delta[index])
        odd[i] = add(odd[i ^ bit], odd_delta[index])
    return tuple(add(x[i], odd[i >> 1] if i & 1 else even[i >> 1])
                 for i in range(16))


def circuit(values, variant='identity'):
    return _circuit(_values(values), lambda a, b: a ^ b,
                    lambda a, b: MUL[a][b], variant)


def circuit_counts(variant='identity'):
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    def mul(a, b):
        counts['fixed_products'] += 1
        return MUL[a][b]
    _circuit((0,) * 16, add, mul, variant)
    return dict(counts)


def _binary_rank(values):
    pivots = [0] * 8
    rank = 0
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if pivots[bit]:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                rank += 1
                break
    return rank


@lru_cache(maxsize=1)
def _half_ranks():
    return tuple(_binary_rank(8 | i for i in range(8) if mask >> i & 1)
                 for mask in range(256))


def split_halves(mask):
    mask = _mask(mask)
    return (sum(((mask >> (2 * i)) & 1) << i for i in range(8)),
            sum(((mask >> (2 * i + 1)) & 1) << i for i in range(8)))


def syndrome_rank(mask):
    """GF2/GF256 rank, using the two affine-three-dimensional halves."""
    e, o = split_halves(mask)
    ranks = _half_ranks()
    return ranks[e] + ranks[o]


def rank_histogram(arity):
    if type(arity) is not int or not 0 <= arity <= 16:
        raise ValueError('arity must be 0..16')
    return dict(sorted(Counter(syndrome_rank(sum(1 << i for i in source))
                               for source in combinations(range(16), arity)).items()))


@lru_cache(maxsize=8)
def _half_quiet_counts(q):
    if type(q) is not int or q < 2 or q & (q - 1):
        raise ValueError('q must be a power of two')
    ranks = _half_ranks()
    values = [q ** (mask.bit_count() - ranks[mask]) for mask in range(256)]
    for bit in range(8):
        for mask in range(256):
            if mask >> bit & 1:
                values[mask] -= values[mask ^ (1 << bit)]
    return tuple(values)


def quiet_support_count(mask, q=256):
    """Exact full-support quiet count on the specified labelled support."""
    e, o = split_halves(mask)
    counts = _half_quiet_counts(q)
    return counts[e] * counts[o]


def quiet_weight_enumerator(q=256):
    """All-amplitude quiet words, by exact total support weight, including0."""
    counts = _half_quiet_counts(q)
    half = [sum(counts[m] for m in range(256) if m.bit_count() == w)
            for w in range(9)]
    return tuple(sum(half[i] * half[w - i]
                     for i in range(max(0, w - 8), min(8, w) + 1))
                 for w in range(17))


def original_quiet_weight_enumerator(q=256):
    """Independent rank-to-weight binomial inversion for the old H*."""
    if type(q) is not int or q < 2 or q & (q - 1):
        raise ValueError('q must be a power of two')
    contained = [0] * 17
    for mask in range(1 << 16):
        w = mask.bit_count()
        contained[w] += q ** (w - original.affine_rank(mask))
    return tuple(sum((-1) ** (w - j) * comb(16 - j, w - j) * contained[j]
                     for j in range(w + 1)) for w in range(17))


@lru_cache(maxsize=240)
def pair_mask_counts(i, j, variant='identity'):
    """Exact255 projective ratio counts; a fixed pair is not resampled."""
    if type(i) is not int or type(j) is not int or not 0 <= i < j < 16:
        raise ValueError('expected ordered distinct source labels')
    h = matrix(variant)
    counts = Counter(sum(bool(h[k][i] ^ MUL[h[k][j]][ratio]) << k
                         for k in range(16)) for ratio in range(1, 256))
    return tuple(sorted(counts.items()))


def pair_summary(variant='identity'):
    classes = Counter()
    worsening = []
    worst = 0
    witness = None
    # Since H is an involution, the same table also classifies the exact
    # input arities that can yield an output of weight2.  None is quiet.
    input_arities_for_weight_two = Counter()
    for i, j in combinations(range(16), 2):
        hist = Counter()
        for mask, count in pair_mask_counts(i, j, variant):
            holes = 16 - mask.bit_count()
            hist[holes] += count
            input_arities_for_weight_two[mask.bit_count()] += count
            if holes > worst:
                worst, witness = holes, (i, j, mask)
        old_worst = max(16 - mask.bit_count()
                        for mask, _ in original.pair_mask_counts(i, j))
        if max(hist) > old_worst:
            worsening.append((i, j, old_worst, max(hist)))
        classes[tuple(sorted(hist.items()))] += 1
    return {
        'ratio_count': 120 * 255,
        'maximum_pair_holes': worst,
        'worst_pair_and_output_mask': witness,
        'pair_law_classes': [{'holes_to_ratio_count': dict(h), 'source_pairs': n}
                             for h, n in sorted(classes.items())],
        'pairs_with_worse_maximum_holes': worsening,
        'input_arities_for_nonquiet_weight_two_output': dict(sorted(input_arities_for_weight_two.items())),
    }


def exact_summary(variant='identity'):
    h, v, b = matrix(variant), syndrome_rows(), coefficient_matrix(variant)
    identity = tuple(tuple(int(i == j) for j in range(16)) for i in range(16))
    gram = matmul(tuple(zip(*v)), v)
    assert all(not x for row in gram for x in row)
    assert determinant(b) != 0
    assert matmul(h, h) == identity
    assert matmul(matmul(v, b), tuple(zip(*v))) == tuple(
        tuple(h[i][j] ^ int(i == j) for j in range(16)) for i in range(16))
    new, old = quiet_weight_enumerator(), original_quiet_weight_enumerator()
    assert sum(new) == 256 ** 8 and sum(old) == 256 ** 11
    assert all(x <= y for x, y in zip(new, old))
    return {
        'construction': construction_identity(variant), 'quadratic_variant': variant,
        'quadratic_coefficient_matrix': quadratic_coefficients(variant),
        'field_order': 256,
        'involution': True, 'zero_entries': sum(not x for row in h for x in row),
        'syndrome_rank': 8, 'coefficient_determinant': determinant(b),
        'rank_histograms': {str(a): rank_histogram(a) for a in (7, 8, 9)},
        'old_rank_histograms': {str(a): original.rank_histogram(a) for a in (7, 8, 9)},
        'quiet_weight_enumerator': new, 'old_quiet_weight_enumerator': old,
        'quiet_minimum_weight': min(w for w in range(1, 17) if new[w]),
        'circuit': circuit_counts(variant), 'old_circuit': original.circuit_counts(),
        'pair_summary': pair_summary(variant),
        'sampler': 'independent uniform nonzero D1/D2 entries, sampled once at setup',
        'whole_code_certificate': False, 'regional_provider_evaluated': False,
        'native_implementation': False, 'old_receipts_transfer': False,
    }


if __name__ == '__main__':
    print(json.dumps(exact_summary(), indent=2))
