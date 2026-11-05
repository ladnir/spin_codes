"""Exact quiet-kernel dimensions for the frozen inverse-column schedule.

Only linear dimensions are counted. No RS incidence probability, global
expurgation, regional operator, or distance certificate is computed here.
"""
from collections import Counter
from itertools import combinations
from math import comb
import time

from hd_bit_syndrome_column_permutations import CONSTRUCTION, PERMUTATIONS
from hd_bit_syndrome_outer import affine_rank


def column_ranks(source_mask):
    if type(source_mask) is not int or not 0 <= source_mask < 65536:
        raise ValueError('A sixteen-coordinate source mask is required')
    rows = tuple(i for i in range(16) if source_mask >> i & 1)
    return tuple(affine_rank(sum(1 << p[i] for i in rows)) for p in PERMUTATIONS)


def dimension_summary(arities=(4, 5, 6, 7, 8, 9), *, max_seconds=10):
    arities = tuple(arities)
    if (not arities or len(set(arities)) != len(arities)
            or any(type(a) is not int or not 1 <= a <= 16 for a in arities)
            or not 0 < max_seconds <= 10):
        raise ValueError('Distinct arities1..16 and a <=10sec budget required')
    started = time.monotonic()
    ranks = tuple(affine_rank(mask) for mask in range(65536))
    images = tuple(tuple(1 << i for i in p) for p in PERMUTATIONS)
    records = []
    for a in arities:
        histogram = Counter()
        hyperplanes = Counter()
        maximum = -1
        witness = None
        processed = 0
        complete = True
        for rows in combinations(range(16), a):
            if processed % 128 == 0 and time.monotonic()-started > max_seconds:
                complete = False
                break
            word_ranks = tuple(ranks[sum(bits[i] for i in rows)] for bits in images)
            dimension = 16*a-sum(word_ranks)
            histogram[dimension] += 1
            if a == 8:
                hyperplanes[word_ranks.count(4)] += 1
            if dimension > maximum:
                maximum = dimension
                witness = dict(source_mask=sum(1 << i for i in rows),
                               source_rows=rows, column_ranks=word_ranks)
            processed += 1
        records.append(dict(arity=a, complete=complete,
            labelled_sets_checked=processed, expected_labelled_sets=comb(16,a),
            quiet_dimension_sum_histogram=dict(sorted(histogram.items())),
            maximum_quiet_dimension_sum=maximum if complete else None,
            rank_four_column_histogram=dict(sorted(hyperplanes.items())) if a == 8 else None,
            maximizing_witness=witness))
        if not complete:
            break
    return dict(construction=CONSTRUCTION, field_order=256, columns=16,
        dimension_formula='sum_j (a-rank_GF256(V_[P_j(A)]))',
        all_requested_arities_complete=len(records)==len(arities) and all(r['complete'] for r in records),
        records=records, elapsed_seconds=time.monotonic()-started,
        source_rows_fixed_across_columns=True, global_setup_probability_computed=False,
        regional_provider_evaluated=False, whole_code_certificate=False)
