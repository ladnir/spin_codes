"""Exact nonlinear input-row schedule; no native implementation or provider.

Column j maps original RS row i to mixer coordinate inv_GF16(i+j), with
inv(0)=0.  GF16 uses the binary label basis and modulus x^4+x+1.
The input RS row set remains fixed across columns.  This schedule is a
new construction, not a permutation symmetry of an existing receipt.
"""
from collections import Counter
from itertools import combinations
import time

from hd_bit_syndrome_outer import affine_rank

CONSTRUCTION = 'G16-Hstar-column-inverse-input-permutations-research'
FIELD_MODULUS = 0x13


def mul(a, b):
    if type(a) is not int or type(b) is not int or not 0 <= a < 16 or not 0 <= b < 16:
        raise ValueError('Binary GF16 labels required')
    answer = 0
    for _ in range(4):
        if b & 1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= FIELD_MODULUS
    return answer


def inverse(a):
    if type(a) is not int or not 0 <= a < 16:
        raise ValueError('Binary GF16 label required')
    if not a:
        return 0
    for b in range(1, 16):
        if mul(a, b) == 1:
            return b
    raise AssertionError('The declared GF16 modulus is not a field')


INVERSE = tuple(inverse(a) for a in range(16))
PERMUTATIONS = tuple(tuple(INVERSE[i ^ j] for i in range(16)) for j in range(16))


def map_mask(source, column):
    if (type(source) is not int or not 0 <= source < 65536
            or type(column) is not int or not 0 <= column < 16):
        raise ValueError('A source mask and column0..15 are required')
    return sum(1 << PERMUTATIONS[column][i] for i in range(16) if source >> i & 1)


def plane_columns(source):
    """Columns permitting a nonzero quiet word on a fixed four-row set."""
    if type(source) is not int or not 0 <= source < 65536 or source.bit_count() != 4:
        raise ValueError('Exactly four fixed original RS rows required')
    indices = tuple(i for i in range(16) if source >> i & 1)
    answer = []
    for j, permutation in enumerate(PERMUTATIONS):
        total = 0
        for i in indices:
            total ^= permutation[i]
        if not total:
            answer.append(j)
    return tuple(answer)


def exact_summary():
    histogram = Counter()
    max_sets = []
    maximum = -1
    for choice in combinations(range(16), 4):
        mask = sum(1 << i for i in choice)
        columns = plane_columns(mask)
        histogram[len(columns)] += 1
        if len(columns) > maximum:
            maximum, max_sets = len(columns), []
        if len(columns) == maximum:
            max_sets.append(dict(source_mask=mask, quiet_columns=columns))
    return dict(construction=CONSTRUCTION, binary_field_modulus=FIELD_MODULUS,
                inverse_table=INVERSE, source_to_mixer_permutations=PERMUTATIONS,
                fixed_source_sets_checked=sum(histogram.values()),
                quiet_column_histogram=dict(sorted(histogram.items())),
                maximum_quiet_columns=maximum,
                minimum_nonquiet_columns_for_four_active_RS_rows=max(0, 9-maximum),
                maximizing_source_sets=max_sets,
                original_RS_length=16, original_RS_minimum_weight=9,
                source_rows_fixed_across_columns=True,
                old_receipts_transfer=False, native_implementation=False,
                regional_provider_evaluated=False, whole_code_certificate=False)


def quiet_activity_incidence(arities=(5, 6, 7), *, max_seconds=55):
    """Exact per-row potential quiet activity, not an RS-word existence test.

    Coordinate i is nonzero in some quiet vector precisely when its column
    is not a coloop.  An all-quiet RS word requires this in at least nine
    columns for EVERY active row.  Failing that necessary condition gives a
    deterministic nonquiet-column floor.  No field amplitudes are sampled.
    """
    arities = tuple(arities)
    if (not arities or any(type(a) is not int or not 1 <= a <= 16 for a in arities)
            or len(set(arities)) != len(arities) or not 0 < max_seconds <= 60):
        raise ValueError('Distinct arities1..16 and a local <=60sec budget required')
    start = time.monotonic()
    ranks = tuple(affine_rank(mask) for mask in range(65536))
    image_bits = tuple(tuple(1 << i for i in row) for row in PERMUTATIONS)
    results = []
    for a in arities:
        histogram = Counter()
        first_by_minimum = {}
        processed = 0
        complete = True
        for choice in combinations(range(16), a):
            if processed % 256 == 0 and time.monotonic() - start > max_seconds:
                complete = False
                break
            activity = [0] * a
            for bits in image_bits:
                mapped = tuple(bits[i] for i in choice)
                mask = sum(mapped)
                rank = ranks[mask]
                for i, bit in enumerate(mapped):
                    activity[i] += ranks[mask ^ bit] == rank
            minimum = min(activity)
            histogram[minimum] += 1
            first_by_minimum.setdefault(minimum, dict(
                source_mask=sum(1 << i for i in choice),
                source_rows=choice, potential_quiet_columns_by_row=activity))
            processed += 1
        results.append(dict(arity=a, complete=complete, labelled_sets_checked=processed,
                            minimum_activity_histogram=dict(sorted(histogram.items())),
                            excluded_all_quiet_sets=sum(n for k, n in histogram.items() if k <= 8),
                            not_excluded_sets=sum(n for k, n in histogram.items() if k >= 9),
                            universal_nonquiet_column_floor=max(0, 9-max(histogram)) if complete else None,
                            first_witness_by_minimum=first_by_minimum))
        if not complete:
            break
    return dict(construction=CONSTRUCTION, elapsed_seconds=time.monotonic()-start,
                all_requested_arities_complete=len(results)==len(arities) and all(r['complete'] for r in results),
                results=results, criterion='every active RS row needs >=9 quiet-capable columns',
                sufficient_for_all_quiet_RS_word=False,
                regional_provider_evaluated=False, whole_code_certificate=False)
