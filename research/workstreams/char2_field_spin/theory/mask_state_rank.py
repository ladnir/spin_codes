"""Bounded actual-state-support rank kernels for packet weights zero to two.

The sixteen states are actual coordinate-support masks, not formal ranks.
Every coefficient is an integer rank/selected-zero count.  Evaluation with
Fraction marks is exact.  Weight three and higher are deliberately unsupported;
there is no unproved distribution of a pooled mass among state masks.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb

import numpy as np

from cauchy_diagonal_packet_reference import LABELS, NODES
from cauchy_joint_return import MUL, INV
from rate_design_coset import power


MULTIPLY = np.asarray(MUL, dtype=np.uint8)
INVERSE = np.asarray(INV, dtype=np.uint8)


def joint_ranks_feasible(equations, following, target_masks,
                         multiply=MULTIPLY, inverse=INVERSE):
    """Exact batched finite-field ranks and linear nonvanishing checks.

    Every current variable is required nonzero.  The following forms must
    have exactly target_masks as their support.  The caller needs q larger
    than the number of required nonzero forms for these separate checks to
    imply simultaneous feasibility.  The GF256 packet has at most ten.
    """
    work = np.array(equations, dtype=np.uint8, copy=True)
    future = np.array(following, dtype=np.uint8, copy=True)
    targets = np.asarray(target_masks)
    batch, rows, columns = work.shape
    indices = np.arange(batch)
    ranks = np.zeros(batch, dtype=np.int64)
    row_indices = np.arange(rows)[None, :]
    for column in range(columns):
        eligible = (row_indices >= ranks[:, None]) & (work[:, :, column] != 0)
        active = eligible.any(axis=1)
        pivots = eligible.argmax(axis=1)
        # Inactive batches do not change, including full-row-rank matrices.
        selected = indices[active]
        if not len(selected):
            continue
        current = ranks[active]
        chosen = pivots[active]
        temporary = work[selected, current].copy()
        work[selected, current] = work[selected, chosen]
        work[selected, chosen] = temporary
        pivot = work[selected, current].copy()
        pivot = multiply[inverse[pivot[:, column]][:, None], pivot]
        work[selected, current] = pivot
        factors = work[selected, :, column].copy()
        factors[np.arange(len(selected)), current] = 0
        work[selected] ^= multiply[factors[:, :, None], pivot[:, None, :]]
        future[selected] ^= multiply[future[selected, :, column, None], pivot[:, None, :]]
        ranks[active] += 1
    # A coordinate variable vanishes on the kernel iff RREF contains its unit row.
    coordinate_alive = ~((work != 0).sum(axis=2) == 1).any(axis=1)
    remaining = (future != 0).any(axis=2)
    actual_masks = (remaining.astype(np.int64) * (1 << np.arange(future.shape[1]))).sum(axis=1)
    return ranks, coordinate_alive & (actual_masks == targets)


def _placement_terms(phase, incoming, inputs):
    """Integer (next mask, rank, selected zeros) counts for fixed placement."""
    inputs = tuple(inputs)
    w = len(inputs)
    if (phase not in range(4) or incoming not in range(16) or w not in (1, 2)
            or tuple(sorted(set(inputs))) != inputs or any(i not in range(16) for i in inputs)):
        raise ValueError('require a valid phase/mask and one or two distinct sorted inputs')
    coordinates = [j for j in range(4) if incoming >> j & 1]
    h = len(coordinates)
    n = h+w
    labels = LABELS[phase]
    output = np.zeros((16, n), dtype=np.uint8)
    for row, u in enumerate(labels):
        for column, j in enumerate(coordinates):
            output[row, column] = MUL[2][power(u ^ NODES[j], 3) ^ 1]
    for column, input_index in enumerate(inputs, h):
        output[input_index, column] = 1
    following = np.zeros((4, n), dtype=np.uint8)
    for column, j in enumerate(coordinates):
        following[j, column] = 1
    for column, input_index in enumerate(inputs, h):
        for j in range(4):
            following[j, column] = power(labels[input_index] ^ NODES[j], 3)

    if not h:
        # Y=X, so selected active positions are impossible and all remaining
        # selected-zero subsets collapse to y**w. Keep their integer counts
        # explicitly so every public profile uses the same mark convention.
        masks = np.arange(16)
        next_forms = np.broadcast_to(following, (16, 4, n))
        zeros = ((masks[:, None] >> np.arange(4)) & 1) == 0
        ranks, alive = joint_ranks_feasible(next_forms*zeros[:, :, None], next_forms, masks)
        return Counter({(int(k), int(r), z): comb(16-w, z)
                        for k, r in zip(masks[alive], ranks[alive]) for z in range(17-w)})

    empty = tuple(j for j in range(16) if j not in inputs)
    active_subsets = [subset for size in range(w+1) for subset in combinations(inputs, size)]
    subsets = []
    for z in range(h):
        for root in combinations(empty, z):
            subsets.extend(root+active for active in active_subsets)
    # Only for ONE input, coordinates outside H cannot cancel each other.
    targets = ([k for k in range(16) if k & (15 ^ incoming) == 15 ^ incoming]
               if w == 1 else list(range(16)))
    zcount, kcount = len(subsets), len(targets)
    maximum_zeros = h-1+w
    base = np.zeros((zcount, maximum_zeros, n), dtype=np.uint8)
    selected_counts = np.empty(zcount, dtype=np.int64)
    for j, subset in enumerate(subsets):
        selected_counts[j] = len(subset)
        if subset:
            base[j, :len(subset)] = output[list(subset)]
    equations = np.zeros((zcount*kcount, maximum_zeros+4, n), dtype=np.uint8)
    equations[:, :maximum_zeros] = np.repeat(base, kcount, axis=0)
    masks = np.tile(np.asarray(targets), zcount)
    next_forms = np.broadcast_to(following, (len(masks), 4, n))
    zero_coordinates = ((masks[:, None] >> np.arange(4)) & 1) == 0
    equations[:, maximum_zeros:] = next_forms*zero_coordinates[:, :, None]
    ranks, alive = joint_ranks_feasible(equations, next_forms, masks)
    zero_sizes = np.repeat(selected_counts, kcount)
    return Counter((int(k), int(r), int(z))
                   for k, r, z in zip(masks[alive], ranks[alive], zero_sizes[alive]))


def _one_placement_terms(phase, incoming, input_index):
    """Compatibility wrapper for the original one-input fixed-placement API."""
    return _placement_terms(phase, incoming, (input_index,))


def translate_mask(mask, offset):
    """GF4 translation permutes the four residue coordinates by XOR."""
    return sum(((mask >> j) & 1) << (j ^ offset) for j in range(4))


def placement_orbits(weight):
    """Representatives and distinct GF4 translates of all input placements."""
    if weight not in (1, 2):
        raise NotImplementedError('only one/two-input placement orbits are supported')
    seen = set()
    result = []
    for inputs in combinations(range(16), weight):
        if inputs in seen:
            continue
        translations = {}
        for offset in range(4):
            target = tuple(sorted(i ^ offset for i in inputs))
            translations.setdefault(target, offset)
        seen.update(translations)
        result.append((inputs, tuple(translations.values())))
    assert len(seen) == comb(16, weight)
    return tuple(result)


@lru_cache(maxsize=12)
def input_profiles(weight, phase=0):
    """Complete 16x16 integer profiles, summing all C(16,weight) placements.

    An entry is a tuple of (rank charge, selected-zero count, integer count).
    Its mark is count*y**(16-z)*(1-y)**z/tau**rank. Divide the total by
    C(16,weight) for the uniform-placement law. Simultaneous GF4 translation
    of labels and residue nodes reduces only exact combinatorial work.
    """
    if phase not in range(4):
        raise ValueError('phase must be zero through three')
    if weight not in (0, 1, 2):
        raise NotImplementedError('only weights zero through two are implemented')
    cells = [[Counter() for _ in range(16)] for _ in range(16)]
    if weight == 0:
        for z in range(17):
            cells[0][0][0, z] = comb(16, z)
        for incoming in range(1, 16):
            for z in range(incoming.bit_count()):
                cells[incoming][incoming][z, z] = comb(16, z)
    else:
        for inputs, offsets in placement_orbits(weight):
            for incoming in range(16):
                terms = _placement_terms(phase, incoming, inputs)
                for offset in offsets:
                    actual_incoming = translate_mask(incoming, offset)
                    for (target, rank, zeros), count in terms.items():
                        actual_target = translate_mask(target, offset)
                        cells[actual_incoming][actual_target][rank, zeros] += count
    return tuple(tuple(tuple((rank, zeros, count) for (rank, zeros), count in sorted(cell.items()))
                       for cell in row) for row in cells)


@lru_cache(maxsize=4)
def one_input_profiles(phase=0):
    """Compatibility profile: original nonzero rows; zero row stays empty."""
    profiles = input_profiles(1, phase)
    return (tuple(() for _ in range(16)),)+profiles[1:]


def placement_profiles(phase, incoming, inputs):
    """Unnormalized integer terms for one fixed input placement and mask."""
    return tuple((target, rank, zeros, count)
                 for (target, rank, zeros), count in sorted(_placement_terms(phase, incoming, inputs).items()))


def local_mask_transfer(weight, y, tau, phase=0):
    """Positive 16-state bound for one packet, conditional on input weight.

    Positive input weights average their fresh uniform input placement. The bound is
    uniform over arbitrary common-history source-form dependencies.  A
    global support sum must still multiply by C(16,weight).
    """
    if weight not in (0, 1, 2):
        raise NotImplementedError('only weights zero through two are implemented')
    if not 0 < y <= 1 or not 1 <= tau <= 255 or phase not in range(4):
        raise ValueError('invalid mark, tilt, or phase')
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    zero = y*0
    result = [[zero for _ in range(16)] for _ in range(16)]
    if weight == 0:
        result[0][0] = y**0
        for incoming in range(1, 16):
            h = incoming.bit_count()
            result[incoming][incoming] = sum(comb(16, z)*y**(16-z)*(1-y)**z/tau**z
                                             for z in range(h))
        return result
    for incoming, row in enumerate(input_profiles(weight, phase)):
        for target, profile in enumerate(row):
            result[incoming][target] = sum(Fraction(count, comb(16, weight))*y**(16-z)*(1-y)**z/tau**rank
                                           for rank, z, count in profile)
    return result


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', type=int, default=0, choices=range(4))
    parser.add_argument('--weight', type=int, default=1, choices=(0, 1, 2))
    args = parser.parse_args()
    value = local_mask_transfer(args.weight, Fraction(9997, 10000), Fraction(35), args.phase)
    print(json.dumps({'phase': args.phase, 'weight': args.weight, 'supported_weights': [0, 1, 2],
                      'nonzero_entries': sum(bool(x) for row in value for x in row),
                      'full_support_return': str(value[15][0]),
                      'whole_code_certificate': False}))
