"""Outward regional placement by truncated matrix-polynomial powering.

Let T_k be a nonnegative local transition envelope for k occupied slots in
a macro containing W slots. A region contains E ordered macros. Define

    B(x) = sum_{k=0}^W C(W,k) T_k x^k.

The average transition envelope for exactly q occupied slots in the region
is [x^q]B(x)^E / C(EW,q). This coefficient sums ordered matrix products.
The matrices T_k need not commute. Binary powering only changes how those
products are parenthesized; it does not exchange any matrix factors.

All powers of x are nonnegative. Discarding degrees above a requested cutoff
therefore cannot affect any retained coefficient, even in later products.
Every coefficient convolution is summed before upward endpoint rounding.
Nonnegative matrices make this rounding monotone under later products.
The final positive binomial division is also rounded upward.

The default backend uses Arb at the caller's current precision. For exact
finite tests, pass fmpq_mat and identity rounding. This helper authenticates
no local model and claims no distance certificate. A proof driver must
authenticate its local envelopes, geometry, source snapshot, and coverage.
"""
from __future__ import annotations

from math import comb

from flint import arb_mat


def rounded(matrix):
    """Replace each finite Arb coefficient by its exact upper endpoint."""
    rows, columns = matrix.nrows(), matrix.ncols()
    values = []
    for i in range(rows):
        row = []
        for j in range(columns):
            value = matrix[i, j]
            if not value.is_finite():
                raise ArithmeticError('finite Arb coefficients required')
            endpoint = value.upper()
            if endpoint < 0:
                raise ArithmeticError('nonnegative coefficient envelopes required')
            row.append(endpoint)
        values.append(row)
    return arb_mat(values)


def _product(left, right, degree, matrix, rounding):
    """Truncated ordered convolution; exactly one callback per completed sum."""
    size = left[0].nrows()
    result = []
    for q in range(min(degree, len(left) + len(right) - 2) + 1):
        total = matrix(size, size)
        first, last = max(0, q - len(right) + 1), min(q, len(left) - 1)
        for j in range(first, last + 1):
            total += left[j] * right[q - j]
        result.append(rounding(total))
    return result


def placement_power(operators, epochs=64, windows=32, matrix=arb_mat,
                    rounding=rounded, maximum_groups=None):
    """Return the 0..D conditional regional operators by binary powering.

    The argument order matches occupancy_model.placement. If maximum_groups
    is omitted, D is the highest available local degree, as in that routine.
    Missing local operators are allowed only above min(D,W). Extra operators
    above W are ignored. Caller-supplied rounding must preserve nonnegative
    upper envelopes; identity rounding is suitable for exact rational input.
    """
    operators = tuple(operators)
    if (type(epochs) is not int or epochs < 1 or type(windows) is not int
            or windows < 1 or not operators or not callable(matrix) or not callable(rounding)):
        raise ValueError('positive integer geometry, local operators, and matrix callbacks required')
    local_degree = min(windows, len(operators) - 1)
    degree = local_degree if maximum_groups is None else maximum_groups
    if (type(degree) is not int or not 0 <= degree <= epochs * windows
            or local_degree < min(degree, windows)):
        raise ValueError('feasible degree and all contributing local occupancy operators required')
    size = operators[0].nrows()
    if not size or any(op.nrows() != size or op.ncols() != size
                       for op in operators[:local_degree + 1]):
        raise ValueError('nonempty square local matrices of the same size required')
    for op in operators[:min(degree, local_degree) + 1]:
        for i in range(size):
            for j in range(size):
                value = op[i, j]
                # Arb inputs must already certify nonnegative local envelopes.
                # The exact-rational backend needs only the same order test.
                if not value >= 0:
                    raise ValueError('local entries must be finite nonnegative envelopes')
                if hasattr(value, 'is_finite') and not value.is_finite():
                    raise ValueError('finite local entries required')
    power = [rounding(operators[k] * comb(windows, k))
             for k in range(min(degree, local_degree) + 1)]
    remaining, result = epochs, None
    while remaining:
        if remaining & 1:
            # Reusing the first selected block avoids multiplying by identity.
            # _product never mutates its input coefficients or lists.
            result = power if result is None else _product(result, power, degree, matrix, rounding)
        remaining >>= 1
        if remaining:
            power = _product(power, power, degree, matrix, rounding)
    return [rounding(value / comb(epochs * windows, q)) for q, value in enumerate(result)]
