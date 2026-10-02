"""Log-domain floating regional proposals, never probability bounds.

All matrix entries are represented by their logarithms. Negative infinity
means an exact structural zero, not a tiny positive value. Positive Arb
upper endpoints are logged before conversion to float. Ordered matrix
products, placement, packet-count mixtures, and final powers stay in the
log domain. No positive entry is removed by an absolute magnitude cutoff.

The conditional placement recurrence appends one macro at a time. For j
occupied slots after e macros, a final occupancy k has probability
C(W,k) C((e-1)W,j-k) / C(eW,j). We retain the ordered product R_(j-k) T_k.
The regional uniform majorant mixes R_j with Bin(q,15/16), then the complete
moment starts at coordinate zero and sums all terminal coordinates. State
is not reset between regions.

Log-sum-exp still uses finite-precision floating arithmetic: negligible
relative contributions can round away when they share a destination.
There are no outward-rounding claims. Every selected witness needs fresh
Arb evaluation with authenticated maps, geometry, and exact outer counts.
This module prepares no maps and writes no files.
"""
from __future__ import annotations

from math import comb, isfinite, log

import numpy as np
from scipy.special import gammaln, logsumexp


def _valid_logs(values, *, family=False):
    values = np.asarray(values, dtype=float)
    dimensions = 3 if family else 2
    if (values.ndim != dimensions or not len(values) or not values.shape[-1]
            or values.shape[-1] != values.shape[-2]
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))):
        raise ValueError('square logarithmic matrices with finite entries or structural -inf required')
    return values


def array_logs(family):
    """Convert finite nonnegative linear arrays; positive subnormals stay positive."""
    values = np.asarray(family, dtype=float)
    if (values.ndim != 3 or not len(values) or not values.shape[1]
            or values.shape[1] != values.shape[2] or not np.isfinite(values).all()
            or np.any(values < 0)):
        raise ValueError('finite nonnegative square operator family required')
    result = np.full(values.shape, -np.inf)
    positive = values > 0
    result[positive] = np.log(values[positive])
    return result


def upper_logs(family):
    """Log positive Arb upper endpoints before converting them to binary64."""
    family = tuple(family)
    if not family or not family[0].nrows():
        raise ValueError('nonempty local Arb family required')
    size = family[0].nrows()
    result = np.full((len(family), size, size), -np.inf)
    for k, matrix in enumerate(family):
        if matrix.nrows() != size or matrix.ncols() != size:
            raise ValueError('consistent square Arb operators required')
        for i in range(size):
            for j in range(size):
                value = matrix[i, j].upper()
                if not value.is_finite() or value < 0:
                    raise ValueError('finite nonnegative local upper endpoints required')
                if value == 0:
                    continue
                logged = float(value.log())
                if not isfinite(logged):
                    raise FloatingPointError('local logarithm exceeds binary64 range')
                result[k, i, j] = logged
    return result


def _multiply(left, right):
    """Ordered log-semiring matrix product; supports a leading batch axis."""
    with np.errstate(over='raise', invalid='raise'):
        terms = left[..., :, :, None] + right[..., None, :, :]
        result = logsumexp(terms, axis=-2)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise FloatingPointError('matrix logarithm exceeded finite range')
    return result


def log_matmul(left, right):
    """Public two-matrix version, preserving the supplied multiplication order."""
    left, right = _valid_logs(left), _valid_logs(right)
    if left.shape != right.shape:
        raise ValueError('matching square logarithmic matrices required')
    return _multiply(left, right)


def log_placement(local_logs, degree, *, epochs, windows=32):
    """Return log R_0 through log R_degree by normalized hypergeometric DP."""
    operators = _valid_logs(local_logs, family=True)
    if (type(degree) is not int or degree < 0 or type(epochs) is not int or epochs < 1
            or type(windows) is not int or windows < 1 or len(operators) != windows+1
            or degree > epochs*windows):
        raise ValueError('complete local family and feasible integer placement geometry required')
    size = operators.shape[1]
    current = np.full((1, size, size), -np.inf)
    np.fill_diagonal(current[0], 0)
    for epoch in range(1, epochs+1):
        limit = min(degree, epoch*windows)
        total, previous = epoch*windows, (epoch-1)*windows
        nxt = np.full((limit+1, size, size), -np.inf)
        for k in range(min(windows, limit)+1):
            js = np.arange(k, min(limit, k+len(current)-1)+1)
            old = js-k
            log_weight = (log(comb(windows, k)) + gammaln(previous+1)
                - gammaln(old+1) - gammaln(previous-old+1) - gammaln(total+1)
                + gammaln(js+1) + gammaln(total-js+1))
            if not np.isfinite(log_weight).all():
                raise FloatingPointError('invalid placement probability logarithm')
            term = _multiply(current[old], operators[k])
            with np.errstate(over='raise', invalid='raise'):
                nxt[js] = np.logaddexp(nxt[js], term + log_weight[:, None, None])
        current = _valid_logs(nxt, family=True)
    return current


def regional_uniform_log(regional_logs, occupancy):
    """Return log E[R_J] for J~Bin(q,15/16), including every feasible J."""
    regional = _valid_logs(regional_logs, family=True)
    if type(occupancy) is not int or not 0 <= occupancy < len(regional):
        raise ValueError('integer occupancy covered by the regional family required')
    q = occupancy
    js = np.arange(q+1)
    weights = (gammaln(q+1)-gammaln(js+1)-gammaln(q-js+1)
               + js*log(15)-q*log(16))
    normalizer = float(logsumexp(weights))
    if not isfinite(normalizer) or abs(normalizer) > 1e-8:
        raise FloatingPointError('binomial logarithms failed floating normalization')
    with np.errstate(over='raise', invalid='raise'):
        result = logsumexp(regional[:q+1] + (weights-normalizer)[:, None, None], axis=0)
    return _valid_logs(result)


def log_power_matrix(log_matrix, exponent):
    """Return log(e_zero M^exponent 1), with continuous state and no reset.

    A structurally zero complete moment is rejected rather than returned as
    an optimistic infinite margin. Zero matrix entries remain permitted.
    """
    matrix = _valid_logs(log_matrix)
    if type(exponent) is not int or exponent < 0:
        raise ValueError('nonnegative integer exponent required')
    if exponent == 0:
        return 0.0
    row = np.full(matrix.shape[0], -np.inf)
    row[0] = 0
    power = matrix
    with np.errstate(over='raise', invalid='raise'):
        while exponent:
            if exponent & 1:
                row = logsumexp(row[:, None]+power, axis=0)
            exponent >>= 1
            if exponent:
                power = _multiply(power, power)
        result = float(logsumexp(row))
    if not isfinite(result):
        raise FloatingPointError('structurally zero or nonfinite complete log moment; use Arb')
    return result
