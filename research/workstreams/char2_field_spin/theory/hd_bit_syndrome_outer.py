"""Exact local algebra for a four-XOR extension of the G16 two-sum mixer.

This is a new outer ensemble, not a transferred distance certificate.
No regional provider, numerical distance bound, or native kernel is run here.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations
import json

from gf8_packet_maps import MUL, determinant, matmul, matvec, xor_sum

CONSTRUCTION = 'G16-RS16x8-D2-bit-syndrome-D1-research'
BASIS = (0x01, 0x0c, 0x50, 0xb0)
OFFSET = 0x06
LABELS = tuple(OFFSET ^ xor_sum(BASIS[j] for j in range(4) if i >> j & 1)
               for i in range(16))


@lru_cache(maxsize=1)
def matrix():
    """H_ij = delta_ij + 1 + u_i u_j + binary dot(i,j), over AES GF256."""
    return tuple(tuple(int(i == j) ^ 1 ^ MUL[LABELS[i]][LABELS[j]]
                       ^ ((i & j).bit_count() & 1)
                       for j in range(16)) for i in range(16))


def affine_matrix(bits=4):
    if not isinstance(bits, int) or not 1 <= bits <= 4:
        raise ValueError('binary label dimension must be 1..4')
    return tuple((1,) + tuple(i >> j & 1 for j in range(bits))
                 for i in range(1 << bits))


def coefficient_matrix():
    """H-I = V B V^T, V=[1,a0,a1,a2,a3], B=I+w w^T."""
    w = (OFFSET,) + BASIS
    return tuple(tuple(int(i == j) ^ MUL[w[i]][w[j]] for j in range(5))
                 for i in range(5))


def _values(values):
    values = tuple(values)
    if len(values) != 16 or any(not isinstance(x, int) or not 0 <= x < 256
                                for x in values):
        raise ValueError('expected sixteen AES-field bytes')
    return values


def direct(values):
    return matvec(matrix(), _values(values))


def _circuit(x, add, mul):
    # Same shared bit-sums as native twoSum16Column; each delta gains A_j.
    p0, p1, p2, p3 = add(x[0], x[1]), add(x[2], x[3]), add(x[4], x[5]), add(x[6], x[7])
    p4, p5, p6, p7 = add(x[8], x[9]), add(x[10], x[11]), add(x[12], x[13]), add(x[14], x[15])
    a0 = add(add(add(x[1], x[3]), add(x[5], x[7])),
             add(add(x[9], x[11]), add(x[13], x[15])))
    a1 = add(add(p1, p3), add(p5, p7))
    q0, q1, q2, q3 = add(p0, p1), add(p2, p3), add(p4, p5), add(p6, p7)
    a2, a3 = add(q1, q3), add(q2, q3)
    s = add(add(q0, q1), a3)
    t = add(a0, add(add(mul(s, OFFSET), mul(a1, BASIS[1])),
                       add(mul(a2, BASIS[2]), mul(a3, BASIS[3]))))
    common = add(s, mul(t, OFFSET))
    delta = (add(t, a0), add(mul(t, BASIS[1]), a1),
             add(mul(t, BASIS[2]), a2), add(mul(t, BASIS[3]), a3))
    corrections = [common] * 16
    for i in range(1, 16):
        bit = i & -i
        corrections[i] = add(corrections[i ^ bit], delta[bit.bit_length() - 1])
    return tuple(add(x[i], corrections[i]) for i in range(16))


def circuit(values):
    return _circuit(_values(values), lambda a, b: a ^ b, lambda a, b: MUL[a][b])


def circuit_counts():
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    def mul(a, b):
        counts['fixed_products'] += 1
        return MUL[a][b]
    _circuit((0,) * 16, add, mul)
    return dict(counts)


def affine_rank(mask, bits=4):
    """Rank of selected augmented binary labels, over GF2 or any extension."""
    if not isinstance(bits, int) or not 1 <= bits <= 4:
        raise ValueError('binary label dimension must be 1..4')
    size = 1 << bits
    if not isinstance(mask, int) or not 0 <= mask < (1 << size):
        raise ValueError('support mask exceeds label set')
    pivots = [0] * (bits + 1)
    answer = 0
    while mask:
        bit = mask & -mask
        value = size | (bit.bit_length() - 1)
        mask ^= bit
        while value:
            lead = value.bit_length() - 1
            if pivots[lead]:
                value ^= pivots[lead]
            else:
                pivots[lead] = value
                answer += 1
                break
    return answer


@lru_cache(maxsize=8)
def quiet_support_counts(q=256, bits=4):
    """Exact forced-NZ counts by support, from binary-matroid Mobius inversion.

    q is a characteristic-two field order. The all-amplitude kernel on S has
    q^(|S|-rank(V_S)); inversion removes every zero coordinate exactly once.
    At bits=4 these are the quiet counts of the proposed H, because B is
    invertible and V has full column rank. Smaller bits are counting toys,
    not a claim that their corresponding update is an involution.
    """
    if not isinstance(q, int) or q < 2 or q & (q - 1):
        raise ValueError('q must be a power of two')
    if not isinstance(bits, int) or not 1 <= bits <= 4:
        raise ValueError('binary label dimension must be 1..4')
    size = 1 << bits
    values = [q ** (mask.bit_count() - affine_rank(mask, bits))
              for mask in range(1 << size)]
    for i in range(size):
        step = 1 << i
        for base in range(0, len(values), 2 * step):
            for j in range(step):
                values[base + step + j] -= values[base + j]
    return tuple(values)


def rank_histogram(arity):
    if not isinstance(arity, int) or not 0 <= arity <= 16:
        raise ValueError('arity must be 0..16')
    return dict(sorted(Counter(affine_rank(sum(1 << i for i in source))
                               for source in combinations(range(16), arity)).items()))


@lru_cache(maxsize=120)
def pair_mask_counts(i, j):
    """Counts of output masks for x_i=1,x_j=ratio in GF256*, not probabilities.

    Fresh nonzero input diagonals make the ratio uniform. The actual pair
    remains fixed across columns. Do not replace these labelled laws by
    independently selected worst pair laws unless an envelope proves it.
    """
    if not isinstance(i, int) or not isinstance(j, int) or not 0 <= i < j < 16:
        raise ValueError('expected ordered distinct source indices 0<=i<j<16')
    h = matrix()
    counts = Counter(sum(bool(h[k][i] ^ MUL[h[k][j]][ratio]) << k
                         for k in range(16)) for ratio in range(1, 256))
    return tuple(sorted(counts.items()))


def exact_summary():
    h, v, b = matrix(), affine_matrix(), coefficient_matrix()
    identity = tuple(tuple(int(i == j) for j in range(16)) for i in range(16))
    gram = matmul(tuple(zip(*v)), v)
    assert all(not x for row in gram for x in row)
    assert determinant(b) != 0
    assert matmul(h, h) == identity
    h_minus_identity = tuple(tuple(h[i][j] ^ int(i == j) for j in range(16))
                             for i in range(16))
    assert matmul(matmul(v, b), tuple(zip(*v))) == h_minus_identity
    classes = Counter()
    worst_holes = 0
    for i, j in combinations(range(16), 2):
        histogram = Counter()
        for mask, count in pair_mask_counts(i, j):
            histogram[16 - mask.bit_count()] += count
        assert sum(histogram.values()) == 255
        worst_holes = max(worst_holes, max(histogram))
        classes[tuple(sorted(histogram.items()))] += 1
    quiet = quiet_support_counts()
    assert min(quiet) >= 0
    assert sum(quiet) == 256 ** 11
    arity8_quiet = Counter(quiet[mask] for mask in range(1 << 16) if mask.bit_count() == 8)
    return {
        'construction': CONSTRUCTION,
        'field_order': 256,
        'labels': LABELS,
        'involution': True,
        'zero_entries': sum(not x for row in h for x in row),
        'syndrome_rank': 5,
        'coefficient_determinant': determinant(b),
        'projective_pairs_checked': 120 * 255,
        'maximum_pair_holes': worst_holes,
        'pair_law_classes': [
            {'holes_to_ratio_count': dict(histogram), 'source_pairs': count}
            for histogram, count in sorted(classes.items())],
        'arity8_rank_histogram': rank_histogram(8),
        'arity8_quiet_count_histogram': dict(sorted(arity8_quiet.items())),
        'quiet_minimum_nonempty_weight': min(mask.bit_count() for mask in range(1, 1 << 16) if quiet[mask]),
        'circuit': circuit_counts(),
        'whole_code_certificate': False,
        'native_implementation': False,
        'distance_provider_evaluated': False,
    }


if __name__ == '__main__':
    print(json.dumps(exact_summary(), indent=2))
