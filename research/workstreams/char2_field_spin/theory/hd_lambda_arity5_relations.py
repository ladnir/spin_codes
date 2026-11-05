"""Exact scheduled Hlambda categories, bounded explicitly through arity five.

Fix an ordered original source tuple A across the sixteen columns. A category
is the GF4-linear relation kernel of its ordered input values. A codimension-d
kernel has prod_{h<d}(q-4**h) raw GF(q) inputs. At q=256, codimension five is
impossible; its zero-count category is retained in the common ordering.

The new arity-five enumeration extends the frozen arity-four provider without
changing its globals. P_j moves original rows before Hlambda; Q_j moves output
rows afterward. Input masks always refer to ORIGINAL source tuple positions.
No common relation kernel across columns, resampling, RS factor, route law,
regional provider, or certificate is asserted here.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product

import hd_quaternary_outer_relations as base
from hd_lambda_outer import matrix as lambda_matrix
from hd_pair_scheduled_outer import SCHEDULES, GF16_POLYNOMIAL, column_permutations

MAX_ARITY = 5
AES_EMBEDDING = base.AES_EMBEDDING
LAMBDA_VALUES = (0xbc, 0xbd)
CONSTRUCTION = 'G16-RS16x8-scheduled-GF4-Hlambda-before-final-diagonal-research'
_AES_DECODE = {value: i for i, value in enumerate(AES_EMBEDDING)}


def _arity(arity):
    if type(arity) is not int or not 0 <= arity <= MAX_ARITY:
        raise ValueError('This isolated interface supports only arity0..5')


def _source(source):
    source = tuple(source)
    if (not 1 <= len(source) <= MAX_ARITY or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('One to five ordered distinct original source labels are required')
    return source


@lru_cache(maxsize=MAX_ARITY+1)
def kernel_subspaces(arity):
    """Canonical RREF bases, ordered by (dimension,packed GF4 row tuple)."""
    _arity(arity)
    if arity <= base.MAX_ARITY:
        return base.kernel_subspaces(arity)
    result = []
    for dimension in range(arity+1):
        for pivots in combinations(range(arity), dimension):
            free = tuple((i, j) for i, pivot in enumerate(pivots)
                         for j in range(pivot+1, arity) if j not in pivots)
            for entries in product(range(4), repeat=len(free)):
                rows = [1 << (2*pivot) for pivot in pivots]
                for (i, j), value in zip(free, entries):
                    rows[i] |= value << (2*j)
                result.append(tuple(rows))
    return tuple(sorted(result, key=lambda rows: (len(rows), rows)))


@lru_cache(maxsize=MAX_ARITY+1)
def kernel_bitsets(arity):
    """Bit v is one exactly when packed GF4 vector v lies in the kernel."""
    _arity(arity)
    if arity <= base.MAX_ARITY:
        return base.kernel_bitsets(arity)
    # Scaling tables avoid repeated digit decoding in the 12,278 kernels.
    scaled = tuple(tuple(sum(base.multiply_gf4(s, row >> (2*i) & 3) << (2*i)
                             for i in range(arity)) for s in range(4))
                   for row in range(4**arity))
    result = []
    for rows in kernel_subspaces(arity):
        members = (0,)
        for row in rows:
            members = tuple(value ^ addition for value in members for addition in scaled[row])
        bits = sum(1 << value for value in members)
        if bits.bit_count() != 4**len(rows):
            raise AssertionError('Canonical RREF enumeration lost linear independence')
        result.append(bits)
    return tuple(result)


@lru_cache(maxsize=32)
def _category_tuples(arity, q):
    _arity(arity)
    base.injection_count(0, q)  # Validate the ambient field without changing base.
    if arity <= base.MAX_ARITY:
        old = base.category_data(arity, q)
        return (old['category_input_masks'], old['category_input_weights'], old['category_counts'])
    masks, weights, counts = [], [], []
    for rows, kernel in zip(kernel_subspaces(arity), kernel_bitsets(arity)):
        mask = sum((not (kernel >> (1 << (2*i)) & 1)) << i for i in range(arity))
        masks.append(mask)
        weights.append(mask.bit_count())
        counts.append(base.injection_count(arity-len(rows), q))
    if sum(counts) != q**arity:
        raise AssertionError('Kernel categories do not partition the ambient input space')
    return tuple(masks), tuple(weights), tuple(counts)


def category_data(arity, q=256):
    """Unmerged exact counts and original-input masks; all categories retained."""
    masks, weights, counts = _category_tuples(arity, q)
    return dict(arity=arity, q=q, coefficient_field_order=4,
                category_order='increasing-dimension-then-canonical-RREF',
                kernel_rrefs=kernel_subspaces(arity), kernel_bitsets=kernel_bitsets(arity),
                category_input_masks=masks, category_input_weights=weights,
                category_counts=counts, raw_allowed_input_mass=q**arity,
                aes_embedding=AES_EMBEDDING if q == 256 else None)


def masks_from_forms(row_forms, arity):
    """Output masks for ordered packed GF4 forms, without merging categories."""
    _arity(arity)
    row_forms = tuple(row_forms)
    if any(type(form) is not int or not 0 <= form < 4**arity for form in row_forms):
        raise ValueError('Every packed GF4 row form must be in range')
    return tuple(sum((not (kernel >> form & 1)) << j for j, form in enumerate(row_forms))
                 for kernel in kernel_bitsets(arity))


def row_forms(source, schedule, column, *, lambda_value=0xbc):
    """Physical output forms of Q_j Hlambda P_j on the original ordered A."""
    source = _source(source)
    if type(lambda_value) is not int or lambda_value not in LAMBDA_VALUES:
        raise ValueError('Use AES lambda0xbc or0xbd, not AES byte2 or3')
    before, after = column_permutations(schedule, column)
    h = lambda_matrix(lambda_value=lambda_value, q=256)
    rows = [0]*16
    for r in range(16):
        rows[after[r]] = sum(_AES_DECODE[h[r][before[i]]] << (2*k)
                             for k, i in enumerate(source))
    return tuple(rows)


def column_masks(source, schedule, column, *, lambda_value=0xbc):
    source = _source(source)
    return masks_from_forms(row_forms(source, schedule, column, lambda_value=lambda_value), len(source))


def column_branches(source, schedule, column, *, lambda_value=0xbc, q=256):
    """Unmerged (original input weight,physical output mask,raw count)."""
    source = _source(source)
    data = category_data(len(source), q)
    masks = column_masks(source, schedule, column, lambda_value=lambda_value)
    return tuple(zip(data['category_input_weights'], masks, data['category_counts']))


def joint_counts(source, schedule, column, *, lambda_value=0xbc, q=256):
    counts = Counter()
    for weight, mask, count in column_branches(source, schedule, column, lambda_value=lambda_value, q=q):
        if count:
            counts[weight, mask] += count
    return dict(sorted(counts.items()))


def column_data(source, schedule, column, *, lambda_value=0xbc, q=256):
    """Single-column data; no source-set, column, or schedule averaging."""
    source = _source(source)
    data = category_data(len(source), q)
    forms = row_forms(source, schedule, column, lambda_value=lambda_value)
    masks = masks_from_forms(forms, len(source))
    before, after = column_permutations(schedule, column)
    return dict(data, original_source=source, original_source_mask=sum(1 << i for i in source),
        mapped_source=tuple(before[i] for i in source), input_permutation=before, output_permutation=after,
        row_forms=forms, category_output_masks=masks,
        branches=tuple(zip(data['category_input_weights'], masks, data['category_counts'])),
        metadata=dict(construction=CONSTRUCTION, lambda_aes_byte=lambda_value,
            lambda_gf4_symbol=AES_EMBEDDING.index(lambda_value), field_order=q,
            coefficient_field_order=4, aes_embedding=AES_EMBEDDING, schedule=schedule,
            column=column, label_field_polynomial=GF16_POLYNOMIAL,
            input_permutation_direction='source i lands at P_j(i) before Hlambda',
            output_permutation_direction='row r lands at Q_j(r) after Hlambda',
            original_A_fixed_across_columns=True, original_source_order_preserved=True,
            categories_unmerged=True, zero_count_categories_retained=True,
            same_relation_kernel_across_columns_required=False, rs_normalization_included=False,
            rowset_resampling=False, schedule_averaging=False, crossed_pair_layer_included=False,
            provider_runs=0, whole_code_certificate=False))
