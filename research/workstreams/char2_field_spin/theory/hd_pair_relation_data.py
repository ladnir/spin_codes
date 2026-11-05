"""Unmerged binary-relation branches for fixed scheduled Hbin source sets.

The original source order identifies the coordinates of F2^a. Every
kernel K of the input map F2^a -> GF(q) is a separate category. Its count
and original-input support do not depend on the column permutations.
Only the output support changes with the prescribed P_j and Q_j.

Source chunks retain their labels across all sixteen columns. Categories
are never merged, including zero-count categories in tiny fields. This
is exact local counting plus FLOAT logarithms, not a regional provider or
distance certificate. The eight crossed-pair blocks are not contracted.
"""
from functools import lru_cache
from itertools import combinations
from math import log2

import numpy as np

from hd_binary_outer_relations import kernel_subspaces, injection_count
from hd_binary_syndrome_outer import matrix
from hd_pair_scheduled_outer import column_permutations


MAX_ARITY = 4
MAX_BRANCH_BATCH = 8192


def _parameters(arity, q):
    if (type(arity) is not int or not 1 <= arity <= MAX_ARITY
            or type(q) is not int or q not in (2, 4, 8, 16, 256)):
        raise ValueError('Arity1..4 and a supported binary-field order required')


@lru_cache(maxsize=4)
def source_sets(arity):
    _parameters(arity, 256)
    return tuple(combinations(range(16), arity))


@lru_cache(maxsize=20)
def category_data(arity, q=256):
    """Stable tuples (kernel_bitset, original_input_mask, weight, count)."""
    _parameters(arity, q)
    result = []
    for kernel in kernel_subspaces(arity):
        dimension = kernel.bit_count().bit_length()-1
        input_mask = sum((not (kernel >> (1 << i) & 1)) << i for i in range(arity))
        count = injection_count(arity-dimension, q)
        result.append((kernel, input_mask, input_mask.bit_count(), count))
    if sum(row[3] for row in result) != q**arity:
        raise AssertionError('Kernel categories do not partition all original inputs')
    return tuple(result)


def _source(source):
    source = tuple(source)
    if (not 1 <= len(source) <= MAX_ARITY or source != tuple(sorted(set(source)))
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('One to four distinct sorted original source rows required')
    return source


def row_forms(source, schedule, column):
    """Actual scheduled binary output forms, indexed by physical output row."""
    source = _source(source)
    before, after = column_permutations(schedule, column)
    h = matrix()
    result = [0]*16
    for row in range(16):
        result[after[row]] = sum(h[row][before[i]] << index for index, i in enumerate(source))
    return tuple(result)


def column_masks(source, schedule, column):
    """One exact support mask per unmerged kernel category, in stable order."""
    source = _source(source)
    forms = np.asarray(row_forms(source, schedule, column), dtype=np.uint16)
    kernels = np.asarray(kernel_subspaces(len(source)), dtype=np.uint16)
    active = 1-((kernels[:, None] >> forms[None, :]) & 1)
    return (active.astype(np.uint32)*(np.uint32(1) << np.arange(16, dtype=np.uint32))).sum(
        axis=1, dtype=np.uint32)


def branch_data(arity, schedule='static', *, sources=None, q=256):
    """Prepare one bounded source chunk in the existing state-kernel layout.

    first/second have shape (8,number_of_A*category_count^2). Shared count,
    weight, and activity arrays have shape (number_of_A,category_count^2).
    At most8192 branch pairs are admitted per block. Use iter_branch_data
    for all a3/a4 row sets. Python-integer counts avoid the q^8 uint64 sum
    overflow at arity4; zero counts have log2 value -inf.
    """
    _parameters(arity, q)
    column_permutations(schedule, 0)
    sources = source_sets(arity) if sources is None else tuple(_source(source) for source in sources)
    if (not sources or any(len(source) != arity for source in sources)
            or len(set(sources)) != len(sources)):
        raise ValueError('A nonempty chunk of distinct original source sets of the selected arity is required')
    categories = category_data(arity, q)
    count = len(categories)
    branches = count*count
    if len(sources)*branches > MAX_BRANCH_BATCH:
        raise ValueError('Source chunk exceeds8192 branch pairs; use iter_branch_data')
    first = np.empty((8, len(sources), branches), dtype=np.uint32)
    second = np.empty_like(first)
    for index, source in enumerate(sources):
        masks = tuple(column_masks(source, schedule, column) for column in range(16))
        for block in range(8):
            first[block, index] = np.repeat(masks[2*block], count)
            second[block, index] = np.tile(masks[2*block+1], count)
    original_masks = np.array([row[1] for row in categories], dtype=np.uint8)
    weights = np.array([row[2] for row in categories], dtype=np.uint8)
    counts = tuple(row[3] for row in categories)
    count_logs = np.array([log2(n) if n else -np.inf for n in counts])
    pair_counts = np.array([a*b for a in counts for b in counts], dtype=object)
    pair_logs = (count_logs[:, None]+count_logs[None, :]).reshape(-1)
    pair_weights = (weights[:, None]+weights[None, :]).reshape(-1)
    pair_activity = ((weights > 0).astype(np.uint8)[:, None]
                     +(weights > 0).astype(np.uint8)[None, :]).reshape(-1)
    pair_union = (original_masks[:, None] | original_masks[None, :]).reshape(-1)
    if sum(pair_counts) != q**(2*arity):
        raise AssertionError('Two-column categories lost raw input mass')
    def shared(values):
        return np.broadcast_to(values, (len(sources), branches))
    return dict(sources=sources, branch_count=branches, category_count=count,
        categories=categories, category_input_masks=original_masks,
        category_input_weights=weights, category_counts=counts,
        first=first.reshape(8, -1), second=second.reshape(8, -1),
        counts=shared(pair_counts), log_counts=shared(pair_logs),
        input_weights=shared(pair_weights), activities=shared(pair_activity),
        input_union_masks=shared(pair_union), deltas=tuple(range(8, 16)),
        columns=tuple((2*j, 2*j+1) for j in range(8)), metadata=metadata(arity, schedule, q=q))


def iter_branch_data(arity, schedule='static', *, batch_size=None, q=256):
    """Yield every original A exactly once, in lexicographic bounded chunks."""
    _parameters(arity, q)
    column_permutations(schedule, 0)
    maximum = MAX_BRANCH_BATCH//len(category_data(arity, q))**2
    batch_size = maximum if batch_size is None else batch_size
    if type(batch_size) is not int or not 1 <= batch_size <= maximum:
        raise ValueError('Batch size must respect the8192-branch limit')
    sources = source_sets(arity)
    for start in range(0, len(sources), batch_size):
        yield branch_data(arity, schedule, sources=sources[start:start+batch_size], q=q)


def metadata(arity, schedule, *, q=256):
    _parameters(arity, q)
    column_permutations(schedule, 0)
    return dict(arity=arity, schedule=schedule, field_order=q,
        categories='all F2-linear relation kernels, unmerged',
        kernel_count=len(category_data(arity, q)), zero_count_categories_retained=True,
        original_source_order_fixed=True, original_A_fixed_across_columns=True,
        input_weight='number of nonzero ORIGINAL source coordinates',
        deterministic_input_and_output_permutations=True, rowset_resampling=False,
        counts_exact_integers=True, logarithms='FLOAT-not-directed',
        regional_provider_runs=0, rs_normalization_included=False, whole_code_certificate=False)
