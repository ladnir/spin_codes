"""Positive FLOAT log-polynomial tail for the paired-column state kernel.

For a fixed labelled row set, branch b has count n_b, original input
weight w_b, and a two-state matrix J_b. Each J_b is the r-th moment of
that branch under the same positive matrix measure. Define, entrywise,

    P_ij(X) = sum_b n_b * J_b[i,j]**(1/r) * X**w_b.

Termwise Holder bounds each coefficient of the original r-group moment
by the corresponding coefficient of P_ij(X)**r. This is an ENTRYWISE
polynomial power, not a matrix power. The resulting block matrices are
then multiplied chronologically. Only afterward do we sum degrees at
least the original-input floor and take the initial-zero row to ones.

The caller retains the same row-set label across blocks. It applies the
fixed-row-set root sum, RS normalization, projective factor, and source
count AFTER this tail. No finite Chernoff input marker or subtraction is
used. All operations are positive ordinary FLOAT in log2, not directed
enclosures. Exact zero is represented by -inf; tiny terms are not floored.
"""
import numpy as np


MAX_DEGREE = 576
MAX_BRANCH_WEIGHT = 4
MAX_GROUPS = 18
MAX_BLOCKS = 8


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite log2 values or exact zero (-inf) required')
    return result


def _polynomial(values):
    result = _logs(values)
    if (result.ndim < 3 or result.shape[-2:] != (2, 2)
            or not 1 <= result.shape[-3] <= MAX_DEGREE+1):
        raise ValueError('Polynomial shape (...,degree,2,2), through degree576, required')
    return result


def branch_polynomial(log_matrices, log_counts, input_weights, groups):
    """Group branches by original input weight before their integer power.

    Inputs have shapes (...,branches,2,2), (...,branches), and
    (...,branches). Output shape is (...,max_weight+1,2,2). Counts may
    include zero (-inf log count); degrees must be integers from0 to4.
    """
    matrices, counts = _logs(log_matrices), _logs(log_counts)
    weights = np.asarray(input_weights)
    if (counts.ndim < 1 or not counts.shape[-1]
            or matrices.shape != counts.shape+(2, 2) or weights.shape != counts.shape
            or not np.issubdtype(weights.dtype, np.integer)
            or np.any(weights < 0) or np.any(weights > MAX_BRANCH_WEIGHT)
            or type(groups) is not int or not 1 <= groups <= MAX_GROUPS):
        raise ValueError('Compatible branches, bounded integer input weights, and groups1..18 required')
    terms = matrices/groups+counts[..., None, None]
    maximum = int(np.max(weights))
    result = np.full(counts.shape[:-1]+(maximum+1, 2, 2), -np.inf)
    for degree in range(maximum+1):
        selected = np.where((weights == degree)[..., None, None], terms, -np.inf)
        result[..., degree, :, :] = np.logaddexp2.reduce(selected, axis=-3)
    return result


def _entrywise_convolution(left, right):
    """Same-entry scalar polynomial products, without summing states."""
    # Internal inputs have already been validated. Iterate over the shorter
    # polynomial; each update is contiguous over degree and both state axes.
    if left.shape[-3] > right.shape[-3]:
        left, right = right, left
    a, b = left.shape[-3], right.shape[-3]
    shape = np.broadcast_shapes(left.shape[:-3], right.shape[:-3])
    result = np.full(shape+(a+b-1, 2, 2), -np.inf)
    for degree in range(a):
        destination = result[..., degree:degree+b, :, :]
        np.logaddexp2(destination,
                      left[..., degree:degree+1, :, :]+right,
                      out=destination)
    return result


def entrywise_power(polynomial, exponent, *, check_time=None):
    """Integer power of EACH matrix entry, preserving every degree."""
    base = _polynomial(polynomial)
    if (type(exponent) is not int or not 0 <= exponent <= MAX_GROUPS
            or exponent*(base.shape[-3]-1) > MAX_DEGREE):
        raise ValueError('Entrywise exponent0..18 and final degree at most576 required')
    # This is entrywise identity1, including the off-diagonal entries.
    answer = np.zeros(base.shape[:-3]+(1, 2, 2))
    # The initial branch polynomial has at most five terms. Multiplying by
    # this fixed short polynomial avoids a costly dense log convolution.
    for _ in range(exponent):
        if check_time is not None:
            check_time()
        answer = _entrywise_convolution(base, answer)
    return answer


def matrix_polynomial_product(left, right, *, check_time=None):
    """Ordered two-state matrix product and positive degree convolution."""
    left, right = _polynomial(left), _polynomial(right)
    a, b = left.shape[-3], right.shape[-3]
    if a+b-2 > MAX_DEGREE:
        raise ValueError('Product degree exceeds576')
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
    """Build the ordered word, retaining the same leading row-set labels.

    blocks is a sequence of at most eight branch-matrix arrays. Shared
    count/weight arrays have the shapes accepted by branch_polynomial.
    """
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
    """Log(e0 * sum_{degree>=minimum} coefficient * ones), per row set."""
    polynomial = _polynomial(polynomial)
    if type(minimum_degree) is not int or not 0 <= minimum_degree <= MAX_DEGREE+1:
        raise ValueError('Nonnegative bounded original-input floor required')
    if minimum_degree >= polynomial.shape[-3]:
        return np.full(polynomial.shape[:-3], -np.inf)
    rows = polynomial[..., minimum_degree:, 0, :]
    return np.logaddexp2.reduce(np.logaddexp2.reduce(rows, axis=-1), axis=-1)


def word_tail(blocks, log_counts, input_weights, groups, minimum_degree, *, check_time=None):
    """Per-fixed-A tail logs; take the A root sum only after this call."""
    word = word_polynomial(blocks, log_counts, input_weights, groups, check_time=check_time)
    return tail_from_polynomial(word, minimum_degree)


def metadata():
    return dict(arithmetic='positive-FLOAT-log2-not-directed',
        degree='original-input-weight', branch_counts_outside_roots=True,
        branch_holder_coefficientwise=True, entrywise_integer_group_power=True,
        chronological_matrix_polynomial_product=True, coefficient_subtraction=False,
        signed_fft=False, finite_input_chernoff_marker=False,
        fixed_rowset_across_blocks=True, rowset_root_sum_included=False,
        rs_normalization_included=False, source_count_included=False,
        original_column_activity_floor_included=False, whole_code_certificate=False)
