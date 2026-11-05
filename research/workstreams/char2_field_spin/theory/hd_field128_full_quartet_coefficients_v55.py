"""Positive full-quartet coefficient calculation; FLOAT proposals only.

Group equal occupancy multisets. Form each chronological matrix product
before its fractional power, then sum all distinct orderings. Histogram
multiplicities remain outside that power. No subtraction or truncation.
"""
from itertools import combinations_with_replacement, islice, permutations
from math import comb, factorial, isfinite, log2

import numpy as np
import hd_field128_gamma_four_packet_positive_v35 as four
import hd_field128_lift_alias_screen as alias


def all_quartet_coefficients(logs, alpha, *, slots=32, chunk=1024,
                             check=lambda: None):
    """Return log coefficients indexed by matrix entry, singles, full pairs.

    A packet occupancy j contributes x^j h_j(y/x^2). Each h_j retains its
    minimum full-count offset. The short positive polynomials let us sum
    65^4 ordered products through C(68,4) unordered occupancy multisets.
    """
    logs = np.asarray(logs, dtype=float)
    if (type(slots) is not int or not 1 <= slots <= 32
            or logs.shape != (2*slots+1, 2, 2)
            or np.any(np.isnan(logs)) or np.any(np.isposinf(logs))
            or not isfinite(alpha) or not 0 < alpha <= 1
            or type(chunk) is not int or not 1 <= chunk <= 4096):
        raise ValueError('Invalid occupancy table, fractional power, or chunk')
    check()
    size, degree = len(logs), 4*slots
    finite = logs[np.isfinite(logs)]
    # This conservative range guard is for ordinary FLOAT coefficients only.
    # A nonzero chronological product has at least one four-factor path.
    if finite.size and (4*alpha*float(finite.min()) < -900
            or 4*alpha*(float(finite.max())+1)+degree*log2(3) > 900):
        raise ValueError('Ordinary positive coefficients need another scale')
    local = []
    for j in range(size):
        low, high = max(0, j-slots), j//2
        values = np.array([comb(slots, j-2*f)*comb(slots-(j-2*f), f)
                           for f in range(low, high+1)], dtype=float)
        local.append((low, values))
    paired = {}
    for a in range(size):
        for b in range(a, size):
            paired[a, b] = (local[a][0]+local[b][0],
                            np.convolve(local[a][1], local[b][1]))
    raw_pairs = four.compose(logs[:, None], logs[None])
    result = np.zeros((4, degree+1, degree+1))
    ordering = tuple(permutations(range(4)))
    stream = combinations_with_replacement(range(size), 4)
    groups = ordered = underflows = 0
    while True:
        check()
        rows = list(islice(stream, chunk))
        if not rows:
            break
        indices = np.asarray(rows, dtype=int)
        weights = np.zeros((len(rows), 2, 2))
        for order in ordering:
            left = raw_pairs[indices[:, order[0]], indices[:, order[1]]]
            right = raw_pairs[indices[:, order[2]], indices[:, order[3]]]
            raw = four.compose(left, right)
            finite_raw = np.isfinite(raw)
            with np.errstate(under='ignore'):
                terms = np.exp2(alpha*raw)
            underflows += int(np.count_nonzero(finite_raw & (terms == 0)))
            if underflows or not np.all(np.isfinite(terms)):
                raise ArithmeticError('Lost a positive raw quartet term')
            weights += terms
        for row, weight in zip(rows, weights):
            divisor = 1
            for j in set(row):
                divisor *= factorial(row.count(j))
            ordered += 24//divisor
            a, b, c, d = row
            low1, first = paired[a, b]
            low2, second = paired[c, d]
            polynomial = np.convolve(first, second)
            f = low1+low2+np.arange(len(polynomial))
            s = sum(row)-2*f
            if np.any(s < 0) or np.any(s+f > degree):
                raise ArithmeticError('Histogram degree/offset mismatch')
            # Within ONE multiset, (s,f) are distinct. No batched repeated
            # advanced-index update is used here.
            result[:, s, f] += (weight.reshape(4)/divisor)[:, None]*polynomial
        groups += len(rows)
    if (groups != comb(size+3, 4) or ordered != size**4
            or not np.all(np.isfinite(result))):
        raise ArithmeticError('Incomplete/nonfinite positive coefficients')
    output = np.full_like(result, -np.inf)
    np.log2(result, out=output, where=result > 0)
    check()
    return output, dict(occupancies=size, multiset_groups=groups,
                       ordered_quartets=ordered, finite_raw_underflows=underflows,
                       maximum_degree=degree, positive_accumulation=True,
                       counts_outside_alpha=True, directed=False)


def point_logs(polynomial, odds):
    degree = np.arange(polynomial.shape[1])[:, None]*odds[0]
    degree = degree+np.arange(polynomial.shape[2])[None, :]*odds[1]
    return np.logaddexp2.reduce((polynomial+degree).reshape((4, -1)),
                                axis=1).reshape((2, 2))


def quartet_alias_table(polynomial, odds, *, grid=128, check=lambda: None):
    """Full nontriangular quartet transform, power32, positive fallback.

    This is not the frozen power64 pair adapter. Degree128 is deliberately
    folded into its residue when grid128 is used; no coefficient is dropped.
    """
    if (polynomial.shape != (4, 129, 129) or grid not in (128, 256)
            or np.any(np.isnan(polynomial)) or np.any(np.isposinf(polynomial))):
        raise ValueError('Expected complete width64 quartet polynomial')
    odds = alias.checked_markers((odds,))[0]
    check()
    point = point_logs(polynomial, odds)
    rho, ratio = alias.perron_logs(point)
    vector = np.array((0., ratio))
    similarity = (vector[None, :]-vector[:, None]).ravel()
    values, lost = alias._frequency(polynomial, odds, similarity-rho, grid)
    rows = values[:, 0, 0].real.reshape((2, 2)).sum(axis=1)
    row_error = float(np.max(abs(rows-1)))
    if not np.all(np.isfinite(values)) or row_error > 1e-9:
        raise ArithmeticError('Normalized quartet row sums failed')
    for _ in range(5):
        check()
        values = alias._matrix_product(values, values)
    residues = np.fft.ifftn(values, axes=(1, 2))
    if not np.all(np.isfinite(residues)):
        raise ArithmeticError('Nonfinite quartet inverse transform')
    imaginary, minimum = float(np.max(abs(residues.imag))), float(np.min(residues.real))
    if imaginary > 1e-10 or minimum < -1e-10:
        raise ArithmeticError('Quartet cancellation diagnostics failed')
    reliable = residues.real >= 1e-10
    logs = np.log2(np.maximum(residues.real, 0)+1e-12)
    logs += 32*rho-similarity[:, None, None]
    positive = four.power(point, 32).ravel()
    check()
    return alias.AliasTable(logs, reliable, positive, dict(
        quartets_per_region=32, normalized_row_sum_error=row_error,
        maximum_imaginary_alias=imaginary, minimum_real_alias=minimum,
        normalized_input_underflow_terms=lost, cancellation_floor=1e-10,
        diagnostic_cushion=1e-12, positive_fallback_required=True,
        directed=False, whole_code_certificate=False))
