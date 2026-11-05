"""Exact scheduled Hlambda columns in stable, unmerged GF4 categories.

Fix an ordered original source tuple A across all sixteen columns. In
column j, source i lands at P_j(i), Hlambda acts, and output row r lands
at Q_j(r). Source order identifies the coordinates of each relation kernel;
neither P_j nor Q_j reorders the original input activity bits.

lambda is the AES byte0xbc or0xbd, not AES byte2. The declared coefficient
field is GF4 with embedding (0,1,0xbc,0xbd). q denotes the ambient field
cardinality for raw input counts. q4 tests use the corresponding symbols2/3.

Counts and original-input masks/weights share one category order across
columns. This does not require an input word to have the same relation
kernel in every column. No categories, row sets, or schedules are averaged.
There is no crossed-pair action, RS factor, regional provider, or bound.
"""
from collections import Counter

from hd_lambda_outer import matrix as lambda_matrix
from hd_pair_scheduled_outer import (
    SCHEDULES, GF16_POLYNOMIAL, column_permutations)
from hd_quaternary_outer_relations import (
    AES_EMBEDDING, category_data as _category_data,
    column_masks as _column_masks, row_forms_from_matrix)

CONSTRUCTION = 'G16-RS16x8-scheduled-GF4-Hlambda-before-final-diagonal-research'
LAMBDA_VALUES = (0xbc, 0xbd)


def _source(source):
    source = tuple(source)
    if (not 1 <= len(source) <= 4 or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('One to four ordered distinct original source labels are required')
    return source


def _lambda(lambda_value):
    if type(lambda_value) is not int or lambda_value not in LAMBDA_VALUES:
        raise ValueError('Use AES lambda0xbc or0xbd, with the explicit GF4 embedding')


def category_data(arity, q=256):
    """Stable raw counts, input masks/weights, kernel RREFs and kernel bitsets.

    The ordering is inherited unchanged from hd_quaternary_outer_relations.
    There are2,7,44,529 categories at arities1,2,3,4. Zero-count categories
    in small ambient fields remain present; no zero-category index is assumed.
    """
    if type(arity) is not int or not 1 <= arity <= 4:
        raise ValueError('The scheduled adapter supports arity1..4')
    return _category_data(arity, q)


def row_forms(source, schedule, column, *, lambda_value=0xbc):
    """Packed GF4 output forms in physical output-row order after Q_j.

    Coefficient k is attached to source[k] BEFORE the input permutation.
    Hlambda is constructed in AES bytes and decoded explicitly to GF4;
    the output row at Q_j(r) receives Hlambda[r,P_j(source[k])].
    """
    source = _source(source)
    _lambda(lambda_value)
    before, after = column_permutations(schedule, column)
    h = lambda_matrix(lambda_value=lambda_value, q=256)
    rows = [None]*16
    for r in range(16):
        rows[after[r]] = tuple(h[r][before[i]] for i in source)
    return row_forms_from_matrix(rows, encoding='aes')


def column_masks(source, schedule, column, *, lambda_value=0xbc):
    """One output mask per stable category, independent of ambient q."""
    source = _source(source)
    forms = row_forms(source, schedule, column, lambda_value=lambda_value)
    return _column_masks(forms, len(source))


def column_branches(source, schedule, column, *, lambda_value=0xbc, q=256):
    """Unmerged (original_input_weight, output_mask, raw_count) tuples.

    This is GF4-kernel order, NOT hd_lambda_outer's hand-written ratio order.
    Equal output masks and zero-count categories are deliberately retained.
    """
    source = _source(source)
    categories = category_data(len(source), q)
    masks = column_masks(source, schedule, column, lambda_value=lambda_value)
    return tuple(zip(categories['category_input_weights'], masks,
                     categories['category_counts']))


def joint_counts(source, schedule, column, *, lambda_value=0xbc, q=256):
    """Optional merged view; stable category data remain available separately."""
    counts = Counter()
    for weight, mask, count in column_branches(source, schedule, column,
                                              lambda_value=lambda_value, q=q):
        if count:
            counts[weight, mask] += count
    return dict(sorted(counts.items()))


def column_data(source, schedule, column, *, lambda_value=0xbc, q=256):
    """Exact single-column data with explicit fixed-source/permutation identity."""
    source = _source(source)
    categories = category_data(len(source), q)
    forms = row_forms(source, schedule, column, lambda_value=lambda_value)
    masks = _column_masks(forms, len(source))
    before, after = column_permutations(schedule, column)
    metadata = dict(construction=CONSTRUCTION, lambda_aes_byte=lambda_value,
        lambda_gf4_symbol=AES_EMBEDDING.index(lambda_value), field_order=q,
        coefficient_field_order=4, aes_embedding=AES_EMBEDDING,
        schedule=schedule, column=column, label_field_polynomial=GF16_POLYNOMIAL,
        input_permutation_direction='source i lands at P_j(i) before Hlambda',
        output_permutation_direction='row r lands at Q_j(r) after Hlambda',
        original_A_fixed_across_columns=True, original_source_order_preserved=True,
        categories_unmerged=True, zero_count_categories_retained=True,
        same_relation_kernel_across_columns_required=False,
        rs_normalization_included=False, rowset_resampling=False,
        schedule_averaging=False, crossed_pair_layer_included=False,
        provider_runs=0, whole_code_certificate=False)
    return dict(categories, original_source=source,
        original_source_mask=sum(1 << i for i in source),
        mapped_source=tuple(before[i] for i in source),
        input_permutation=before, output_permutation=after,
        row_forms=forms, category_output_masks=masks,
        branches=tuple(zip(categories['category_input_weights'], masks,
                           categories['category_counts'])), metadata=metadata)
