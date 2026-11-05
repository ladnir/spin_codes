"""Exact inside cosets conditioned on one nonzero outside output.

This is local counting for the fixed AES-field bit-syndrome matrix H*.
It is not a regional operator or a distance certificate.  A is one fixed
labelled eight-row input set; j is outside A.  Counts are for y_j=1 and
y_i=0 at the other outside positions.  Every other nonzero value at j
has the same support counts, so summing those values multiplies by 255.

Write H*=I+V B V^T, with V_i=(1, bits(i)) and B=I+w w^T.  A feasible
outside singleton determines a binary c from V_outside c=e_j.  Then
V_A^T y_A=d=B^-1 c+V_A^T V_A c=B^-1 c+V_j.  The last equality uses
V^T V=0.  For S contained in A, the number supported within S is
256^(|S|-rank V_S) if d belongs to span_F256(V_S), and zero otherwise.
Subset inversion gives the exact support counts.

Binary affine symmetries are NOT assumed to preserve H*.  All 240
possible (c,j) syndromes, and all 1920 labelled A, are checked explicitly.
"""
from collections import Counter
from functools import lru_cache
import json

from gf8_packet_maps import INV, MUL, binary_rank, xor_sum
from hd_bit_syndrome_outer import (
    BASIS, CONSTRUCTION, OFFSET, affine_matrix, affine_rank,
)

Q = 256
V = affine_matrix()
W = (OFFSET,) + BASIS
DELTA = 1 ^ xor_sum(MUL[x][x] for x in W)
GENERIC_WEIGHTS = (0, 0, 0, 0, 0, 28, 5292, 448084, 16323812)
EXCEPTION_WEIGHTS = (0, 0, 0, 0, 4, 268, 4548, 448836, 16323560)
GENERIC_CAPS = (0, 0, 0, 0, 0, 1, 252, 64012, 16323812)
EXCEPTION_CAPS = (0, 0, 0, 0, 1, 252, 253, 64264, 16323560)
UNIVERSAL_CAPS = tuple(map(max, zip(GENERIC_CAPS, EXCEPTION_CAPS)))


def _source(source, outside_position):
    source = tuple(source)
    if (len(source) != 8 or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 16 for i in source)
            or type(outside_position) is not int
            or not 0 <= outside_position < 16 or outside_position in source):
        raise ValueError('Sorted labelled A of size8 and j outside A required')
    return source


def _dot_binary(mask, values):
    return xor_sum(value for bit, value in enumerate(values) if mask >> bit & 1)


@lru_cache(maxsize=65536)
def _rank(mask):
    return affine_rank(mask)


def syndrome(c_mask, outside_position):
    """Return the exact five-byte syndrome for a normalized singleton."""
    if (type(c_mask) is not int or not 2 <= c_mask < 32
            or type(outside_position) is not int or not 0 <= outside_position < 16
            or _dot_binary(c_mask, V[outside_position]) != 1):
        raise ValueError('Nonconstant binary affine form with f(j)=1 required')
    c = tuple(c_mask >> bit & 1 for bit in range(5))
    beta = MUL[_dot_binary(c_mask, W)][INV[DELTA]]
    d = tuple(a ^ b ^ MUL[beta][w] for a, b, w in zip(c, V[outside_position], W))
    r = binary_rank(d)
    if r not in (4, 5):
        raise AssertionError('The fixed byte-matrix syndrome rank changed')
    annihilator = next((h for h in range(1, 32) if not _dot_binary(h, d)), 0)
    if (r == 5) != (annihilator == 0):
        raise AssertionError('Syndrome annihilator disagrees with binary rank')
    plane = (sum((not _dot_binary(annihilator, v)) << i for i, v in enumerate(V))
             if annihilator else 0)
    return dict(c_mask=c_mask, c=c, outside_position=outside_position,
                inside_syndrome=d, binary_rank=r, annihilator=annihilator,
                annihilator_position_mask=plane)


def singleton_data(source, outside_position):
    """Return coset data, or None if this outside singleton is impossible."""
    source = _source(source, outside_position)
    outside = tuple(i for i in range(16) if i not in source)
    matches = [c for c in range(2, 32)
               if all(_dot_binary(c, V[i]) == int(i == outside_position) for i in outside)]
    if not matches:
        return None
    if len(matches) != 1 or _rank(sum(1 << i for i in source)) != 5:
        raise AssertionError('Feasible singleton must determine one rank5 coset')
    return dict(source=source, outside=outside, **syndrome(matches[0], outside_position))


def _within_counts(source, data):
    masks = [0] * 256
    counts = [0] * 256
    plane = data['annihilator_position_mask']
    for local in range(1, 256):
        bit = local & -local
        masks[local] = masks[local ^ bit] | (1 << source[bit.bit_length() - 1])
        r = _rank(masks[local])
        # Binary column spans of dimension5 contain every target.  When d
        # has binary rank4 its bit planes span exactly ker(annihilator).
        if r == 5 or (r == 4 and plane and not (masks[local] & ~plane)):
            counts[local] = Q ** (local.bit_count() - r)
    return masks, counts


def counts_for(source, outside_position):
    """Exact {GLOBAL inside-support mask:count}, normalized outside value1."""
    data = singleton_data(source, outside_position)
    if data is None:
        return {}
    masks, counts = _within_counts(data['source'], data)
    for bit in range(8):
        step = 1 << bit
        for base in range(0, 256, 2 * step):
            for offset in range(step):
                counts[base + step + offset] -= counts[base + offset]
    if min(counts) < 0 or sum(counts) != Q ** 3:
        raise AssertionError('Exact affine-fiber count identity failed')
    return {mask: count for mask, count in zip(masks, counts) if count}


def weight_counts(counts):
    values = [0] * 9
    for mask, count in counts.items():
        values[mask.bit_count()] += count
    return tuple(values)


def support_caps(counts):
    values = [0] * 9
    for mask, count in counts.items():
        w = mask.bit_count()
        values[w] = max(values[w], count)
    return tuple(values)


def all_singletons():
    """Yield each of the1920 feasible labelled (A,j) exactly once."""
    for c in range(2, 32):
        values = tuple(_dot_binary(c, v) for v in V)
        for j in range(16):
            if not values[j]:
                continue
            for k in range(16):
                if not values[k]:
                    yield tuple(i for i in range(16) if (values[i] and i != j) or i == k), j


def audit_all():
    """Reproduce all exact local counts; no regional provider or code gate."""
    ranks = Counter(syndrome(c, j)['binary_rank'] for c in range(2, 32)
                    for j in range(16) if _dot_binary(c, V[j]))
    weight_hist, cap_hist = Counter(), Counter()
    seen, exceptions = set(), []
    for source, j in all_singletons():
        mask = sum(1 << i for i in source)
        if mask in seen:
            raise AssertionError('A feasible labelled row set was duplicated')
        seen.add(mask)
        counts = counts_for(source, j)
        weights, caps = weight_counts(counts), support_caps(counts)
        weight_hist[weights] += 1
        cap_hist[caps] += 1
        if weights[4]:
            exceptions.append(dict(source=source, outside_position=j,
                                   **{k: v for k, v in singleton_data(source, j).items()
                                      if k in ('c', 'inside_syndrome', 'annihilator')}))
    if (len(seen) != 1920 or ranks != {5: 224, 4: 16}
            or weight_hist != {GENERIC_WEIGHTS: 1888, EXCEPTION_WEIGHTS: 32}
            or cap_hist != {GENERIC_CAPS: 1888, EXCEPTION_CAPS: 32}):
        raise AssertionError('The fixed-matrix finite audit changed')
    return dict(construction=CONSTRUCTION, field_order=Q, outside_value=1,
                outside_nonzero_value_multiplier=255, rowsets=1920,
                syndrome_cases=240, syndrome_binary_ranks=dict(ranks),
                generic_fibers=1888, exceptional_fibers=32,
                generic_weight_counts=GENERIC_WEIGHTS,
                exception_weight_counts=EXCEPTION_WEIGHTS,
                generic_per_support_caps=GENERIC_CAPS,
                exception_per_support_caps=EXCEPTION_CAPS,
                universal_per_support_caps=UNIVERSAL_CAPS,
                exceptions=exceptions, total_output_weight_minimum=5,
                regional_provider=False, whole_code_certificate=False)


if __name__ == '__main__':
    print(json.dumps(audit_all(), indent=2))
