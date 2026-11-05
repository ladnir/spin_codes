"""Exact low-arity branch data for deterministic Hbin row schedules.

Fix an original RS row set A for all sixteen columns. In column j, source
row i is sent to P_j(i), Hbin is applied, and output row r is sent to Q_j(r).
The crossed-pair layer follows these permutations; it is not included here.
Input weight always counts the original pre-mixer symbols.

This module prepares exact counts and labelled support masks for arities
one and two. It does not resample A, average schedules, run a regional
provider, or compute a distance bound. Log counts are ordinary FLOAT; the
integer counts are also returned. Existing state-kernel/tail APIs accept
each block's flattened first/second arrays and the shared count arrays.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import log2

import numpy as np

from hd_binary_syndrome_outer import column_mask, two_row_joint_counts


SCHEDULES = ('static', 'input-xor', 'input-inverse', 'output-xor',
             'input-inverse-output-xor')
GF16_POLYNOMIAL = 0x13


def gf16_multiply(left, right):
    """GF16 labels in the polynomial basis modulo x^4+x+1."""
    if type(left) is not int or type(right) is not int or not 0 <= left < 16 or not 0 <= right < 16:
        raise ValueError('GF16 labels0..15 required')
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 16:
            left ^= GF16_POLYNOMIAL
    return result


@lru_cache(maxsize=1)
def inverse_table():
    """The prescribed permutation inv(0)=0, inv(x)=x^14 otherwise."""
    result = [0]
    for value in range(1, 16):
        power = 1
        for _ in range(14):
            power = gf16_multiply(power, value)
        if gf16_multiply(value, power) != 1:
            raise AssertionError('GF16 inverse identity failed')
        result.append(power)
    return tuple(result)


@lru_cache(maxsize=80)
def column_permutations(schedule, column):
    """Return source-to-destination permutations (P_j,Q_j)."""
    if schedule not in SCHEDULES or type(column) is not int or not 0 <= column < 16:
        raise ValueError('A prescribed schedule and column0..15 are required')
    identity = tuple(range(16))
    if schedule == 'input-xor':
        before = tuple(i ^ column for i in range(16))
    elif schedule in ('input-inverse', 'input-inverse-output-xor'):
        inverses = inverse_table()
        before = tuple(inverses[i ^ column] for i in range(16))
    else:
        before = identity
    after = (tuple(i ^ column for i in range(16))
             if schedule in ('output-xor', 'input-inverse-output-xor') else identity)
    return before, after


def permute_mask(mask, permutation):
    """Send each set bit i to permutation[i], without changing its weight."""
    permutation = tuple(permutation)
    if (type(mask) is not int or not 0 <= mask < 65536
            or len(permutation) != 16 or set(permutation) != set(range(16))):
        raise ValueError('A 16-bit support mask and row permutation required')
    result = 0
    while mask:
        low = mask & -mask
        result |= 1 << permutation[low.bit_length()-1]
        mask ^= low
    return result


def column_branches(source, schedule, column, *, q=256):
    """Return ordered (original input weight, output mask, exact count).

    The order is zero, singleton original-row0, singleton original-row1,
    equal nonzero inputs, distinct nonzero inputs. A singleton source uses
    only the first two categories. Input permutations never reorder these
    categories, even when P_j reverses the two source labels.
    """
    source = tuple(source)
    if (len(source) not in (1, 2) or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 16 for i in source)
            or type(q) is not int or q not in (4, 8, 16, 256)):
        raise ValueError('One/two sorted original rows and supported binary field order required')
    before, after = column_permutations(schedule, column)
    mapped = tuple(before[i] for i in source)
    masks = tuple(permute_mask(column_mask(i), after) for i in mapped)
    branches = [(0, 0, 1), (1, masks[0], q-1)]
    if len(source) == 2:
        branches.extend(((1, masks[1], q-1), (2, masks[0] ^ masks[1], q-1),
                         (2, masks[0] | masks[1], (q-1)*(q-2))))
        # Compare the stable category representation with the established
        # exact joint-count helper, after the actual output permutation.
        exact = two_row_joint_counts(*sorted(mapped), q=q)
        expected = Counter()
        for (weight, mask), count in exact.items():
            expected[weight, permute_mask(mask, after)] += count
        observed = Counter()
        for weight, mask, count in branches:
            observed[weight, mask] += count
        if observed != expected:
            raise AssertionError('Scheduled branches disagree with the exact joint-count helper')
    if sum(count for _, _, count in branches) != q**len(source):
        raise AssertionError('Original-input mass changed')
    return tuple(branches)


def branch_data(arity, schedule='static', *, q=256):
    """Prepare all original A, with per-block supports and invariant counts.

    first/second have shape (8, number_of_A*branch_count). Their block b
    uses columns2b and2b+1. Counts, weights, and activities have shape
    (number_of_A,branch_count), identical across the eight blocks.
    """
    if type(arity) is not int or arity not in (1, 2):
        raise ValueError('Only arities1 and2 are prepared')
    column_permutations(schedule, 0)
    if type(q) is not int or q not in (4, 8, 16, 256):
        raise ValueError('Supported binary-field order required')
    sources = tuple(combinations(range(16), arity))
    categories = 2 if arity == 1 else 5
    branch_count = categories**2
    first = np.empty((8, len(sources), branch_count), dtype=np.uint32)
    second = np.empty_like(first)
    counts = np.empty((len(sources), branch_count), dtype=np.uint64)
    weights = np.empty((len(sources), branch_count), dtype=np.uint8)
    activities = np.empty_like(weights)
    # The logarithm follows the existing static caller's sum-of-two-logs
    # convention; integer counts remain available independently.
    log_counts = np.empty((len(sources), branch_count))
    for source_index, source in enumerate(sources):
        columns = tuple(column_branches(source, schedule, j, q=q) for j in range(16))
        for block in range(8):
            for index, (left, right) in enumerate(product(columns[2*block], columns[2*block+1])):
                first[block, source_index, index] = left[1]
                second[block, source_index, index] = right[1]
                count, weight = left[2]*right[2], left[0]+right[0]
                activity = int(left[0] > 0)+int(right[0] > 0)
                if block == 0:
                    counts[source_index, index] = count
                    weights[source_index, index] = weight
                    activities[source_index, index] = activity
                    log_counts[source_index, index] = log2(left[2])+log2(right[2])
                elif (count != counts[source_index, index] or weight != weights[source_index, index]
                      or activity != activities[source_index, index]):
                    raise AssertionError('Permutation changed a branch category')
    return dict(sources=sources, branch_count=branch_count,
        first=first.reshape(8, -1), second=second.reshape(8, -1),
        counts=counts, log_counts=log_counts, input_weights=weights, activities=activities,
        columns=tuple((2*b, 2*b+1) for b in range(8)), deltas=tuple(range(8, 16)),
        input_permutations=tuple(column_permutations(schedule, j)[0] for j in range(16)),
        output_permutations=tuple(column_permutations(schedule, j)[1] for j in range(16)),
        metadata=metadata(schedule, q=q))


def metadata(schedule, *, q=256):
    column_permutations(schedule, 0)
    if type(q) is not int or q not in (4, 8, 16, 256):
        raise ValueError('Supported binary-field order required')
    return dict(schedule=schedule, field_order=q, mixer='Hbin=I+J+four-binary-syndromes',
        permutation_field='GF16-polynomial-0x13', inverse_zero=0,
        input_permutation_direction='source i lands at P_j(i) BEFORE Hbin',
        output_permutation_direction='Hbin row r lands at Q_j(r) BEFORE crossed pairs',
        original_A_fixed_across_all_columns=True, deterministic_schedule=True,
        rowset_resampling=False, schedule_averaging=False, original_input_weight_retained=True,
        counts_exact_integers=True, logarithms='FLOAT-not-directed',
        full_independent_D1_comparison=True, rs_normalization_included=False,
        provider_runs=0, native_implementation_changed=False, whole_code_certificate=False)
