"""Bounded positive log-polynomial tails for binary outer arities1..4.

This extends the reviewed hd_pair_state_tail interface to branch input
weight8 and word degree1152. It does not alter that frozen module. The
coefficientwise Holder argument, entrywise group powers, chronological
matrix products, and final fixed-A tail are unchanged.

Only bound-independent log validation and scalar positive convolution
are reused from the frozen source. No module globals are patched. The
caller applies the row-set root sum and all external normalizers.
Ordinary FLOAT log2 arithmetic is not a directed certificate.
"""
import numpy as np

from hd_pair_state_tail import _logs, _entrywise_convolution


MAX_DEGREE = 1152
MAX_BRANCH_WEIGHT = 8
MAX_GROUPS = 18
MAX_BLOCKS = 8


def _polynomial(values):
    result = _logs(values)
    if (result.ndim < 3 or result.shape[-2:] != (2, 2)
            or not 1 <= result.shape[-3] <= MAX_DEGREE+1):
        raise ValueError('Polynomial shape (...,degree,2,2), through degree1152, required')
    return result


def branch_polynomial(log_matrices, log_counts, input_weights, groups):
    """Coefficients of sum_b n_b J_b**(1/r) X**original_weight_b."""
    matrices, counts = _logs(log_matrices), _logs(log_counts)
    weights = np.asarray(input_weights)
    if (counts.ndim < 1 or not counts.shape[-1]
            or matrices.shape != counts.shape+(2, 2) or weights.shape != counts.shape
            or not np.issubdtype(weights.dtype, np.integer)
            or np.any(weights < 0) or np.any(weights > MAX_BRANCH_WEIGHT)
            or type(groups) is not int or not 1 <= groups <= MAX_GROUPS):
        raise ValueError('Compatible branches, integer input weights0..8, and groups1..18 required')
    terms = matrices/groups+counts[..., None, None]
    maximum = int(np.max(weights))
    result = np.full(counts.shape[:-1]+(maximum+1, 2, 2), -np.inf)
    for degree in range(maximum+1):
        selected = np.where((weights == degree)[..., None, None], terms, -np.inf)
        result[..., degree, :, :] = np.logaddexp2.reduce(selected, axis=-3)
    return result


def entrywise_power(polynomial, exponent, *, check_time=None):
    """Integer positive polynomial power of each matrix entry separately."""
    base = _polynomial(polynomial)
    if (type(exponent) is not int or not 0 <= exponent <= MAX_GROUPS
            or exponent*(base.shape[-3]-1) > MAX_DEGREE):
        raise ValueError('Entrywise exponent0..18 and final degree at most1152 required')
    answer = np.zeros(base.shape[:-3]+(1, 2, 2))
    for _ in range(exponent):
        if check_time is not None:
            check_time()
        answer = _entrywise_convolution(base, answer)
    return answer


def matrix_polynomial_product(left, right, *, check_time=None):
    """Ordered two-state product with positive log-degree convolution."""
    left, right = _polynomial(left), _polynomial(right)
    a, b = left.shape[-3], right.shape[-3]
    if a+b-2 > MAX_DEGREE:
        raise ValueError('Product degree exceeds1152')
    shape = np.broadcast_shapes(left.shape[:-3], right.shape[:-3])
    result = np.full(shape+(a+b-1, 2, 2), -np.inf)
    if a <= b:
        for degree in range(a):
            if check_time is not None:
                check_time()
            fixed = left[..., degree:degree+1, :, :]
            product = np.logaddexp2(
                fixed[..., :, 0, None]+right[..., None, 0, :],
                fixed[..., :, 1, None]+right[..., None, 1, :])
            destination = result[..., degree:degree+b, :, :]
            np.logaddexp2(destination, product, out=destination)
    else:
        for degree in range(b):
            if check_time is not None:
                check_time()
            fixed = right[..., degree:degree+1, :, :]
            product = np.logaddexp2(
                left[..., :, 0, None]+fixed[..., None, 0, :],
                left[..., :, 1, None]+fixed[..., None, 1, :])
            destination = result[..., degree:degree+a, :, :]
            np.logaddexp2(destination, product, out=destination)
    return _logs(result)


def word_polynomial(blocks, log_counts, input_weights, groups, *, check_time=None):
    """Same fixed-A labels across one to eight chronological blocks."""
    if not hasattr(blocks, '__len__') or not 1 <= len(blocks) <= MAX_BLOCKS:
        raise ValueError('One to eight chronological blocks required')
    counts = _logs(log_counts)
    if counts.ndim < 1:
        raise ValueError('A branch axis is required')
    word = np.full(counts.shape[:-1]+(1, 2, 2), -np.inf)
    word[..., 0, 0, 0] = word[..., 0, 1, 1] = 0.
    for block in blocks:
        if check_time is not None:
            check_time()
        polynomial = branch_polynomial(block, counts, input_weights, groups)
        current = entrywise_power(polynomial, groups, check_time=check_time)
        word = matrix_polynomial_product(word, current, check_time=check_time)
    return word


def tail_from_polynomial(polynomial, minimum_degree):
    """Log(e0 * sum_{degree>=minimum} coefficient * ones), per fixed A."""
    polynomial = _polynomial(polynomial)
    if type(minimum_degree) is not int or not 0 <= minimum_degree <= MAX_DEGREE+1:
        raise ValueError('Nonnegative bounded original-input floor required')
    if minimum_degree >= polynomial.shape[-3]:
        return np.full(polynomial.shape[:-3], -np.inf)
    rows = polynomial[..., minimum_degree:, 0, :]
    return np.logaddexp2.reduce(np.logaddexp2.reduce(rows, axis=-1), axis=-1)


def word_tail(blocks, log_counts, input_weights, groups, minimum_degree, *, check_time=None):
    word = word_polynomial(blocks, log_counts, input_weights, groups, check_time=check_time)
    return tail_from_polynomial(word, minimum_degree)


def metadata():
    return dict(arithmetic='positive-FLOAT-log2-not-directed',
        degree='original-input-weight', maximum_branch_weight=8, maximum_word_degree=1152,
        branch_counts_outside_roots=True, branch_holder_coefficientwise=True,
        entrywise_integer_group_power=True, chronological_matrix_polynomial_product=True,
        coefficient_subtraction=False, signed_fft=False, fixed_rowset_across_blocks=True,
        rowset_root_sum_included=False, rs_normalization_included=False,
        source_count_included=False, original_column_activity_floor_included=False,
        whole_code_certificate=False)
