"""Exact joint pre-mixer weight/output-mask counts for the affine two-sum H.

This is a counting helper, not a distance certificate.  It strengthens the
relaxation of a fixed set of a nonzero RS[16,8] rows: their total pre-H
support across columns is at least9a, not merely9 nonempty columns.
"""
from functools import lru_cache
from fractions import Fraction as F
from math import log2

import numpy as np

from hd_two_sum_closed_kernel import branch_counts, choose


@lru_cache(None)
def ambient_count(arity, overlap, weight, group=16):
    """N_A(T) from the sizes a, |A intersect T| and |T|."""
    if (group not in (8, 16) or any(type(x) is not int for x in (arity, overlap, weight))
            or not 0 <= arity <= group or not 0 <= weight <= group
            or not max(0, arity+weight-group) <= overlap <= min(arity, weight)):
        raise ValueError('Feasible fixed input/output supports required')
    quiet, hole, full = branch_counts(arity, overlap, group)
    outside = group-arity
    selected = weight-overlap
    return ((quiet if selected == 0 else 0)
            +(hole if outside-selected == 1 else 0)
            +(full if selected == outside else 0))


@lru_cache(None)
def exact_support_count(arity, overlap, weight, group=16):
    """Count inputs with support EXACTLY A and outputs with support T."""
    ambient_count(arity, overlap, weight, group)  # Validate before summing.
    result = sum((-1)**(arity-inside-outside)*choose(overlap, inside)
                 *choose(arity-overlap, outside)
                 *ambient_count(inside+outside, inside, weight, group)
                 for inside in range(overlap+1)
                 for outside in range(arity-overlap+1))
    if result < 0:
        raise ArithmeticError('Negative exact joint-support count')
    return result


@lru_cache(None)
def input_weight_counts(arity, overlap, weight, group=16):
    """Positive coefficients Q[k] of sum_x v^wt(x), with supp(x)<=A."""
    ambient = ambient_count(arity, overlap, weight, group)
    result = tuple(sum(choose(overlap, inside)*choose(arity-overlap, k-inside)
                       *exact_support_count(k, inside, weight, group)
                       for inside in range(max(0, k-(arity-overlap)), min(k, overlap)+1))
                   for k in range(arity+1))
    if sum(result) != ambient or min(result) < 0:
        raise ArithmeticError('Joint-weight polynomial failed its unmarked identity')
    return result


def marked_count(arity, overlap, weight, marker=F(1), group=16):
    marker = F(marker)
    if marker < 1:
        raise ValueError('Pre-H lower-support marker must be at least1')
    return sum((coefficient*marker**k for k, coefficient in enumerate(
        input_weight_counts(arity, overlap, weight, group))), F(0))


def type_weights(alpha, log_marker, *, group=16, columns=16, dimension=8, minimum=9):
    """FLOAT W_a(T), including v^(-minimum*a*alpha/columns).

    The support marker is introduced BEFORE fractionalization, so it IS
    alpha-powered. Row choices remain outside alpha inside the Holder sum.
    One marker is fixed for each type across all paths and columns. There
    is no separate nonempty-column marker in this helper.
    """
    alpha = F(alpha)
    if (not 0 < alpha <= 1 or not np.isfinite(log_marker) or log_marker < 0
            or group not in (8, 16) or type(columns) is not int or columns < 1
            or not 0 < dimension < columns or not 1 <= minimum <= columns):
        raise ValueError('Valid fractional exponent and pre-H support marker required')
    result = {}
    for arity in range(1, group+1):
        weights = []
        for weight in range(group+1):
            terms = []
            for overlap in range(max(0, arity+weight-group), min(arity, weight)+1):
                coefficients = input_weight_counts(arity, overlap, weight, group)
                marked = [log2(c)+k*log_marker for k, c in enumerate(coefficients) if c]
                if marked:
                    mass = float(np.logaddexp2.reduce(marked))
                    multiplicity = choose(weight, overlap)*choose(group-weight, arity-overlap)
                    terms.append(log2(multiplicity)+columns*float(alpha)*mass)
            total = float(np.logaddexp2.reduce(terms)) if terms else -np.inf
            weights.append((total-(columns-dimension)*arity*float(alpha)*log2(255)
                            -minimum*arity*float(alpha)*log_marker)/columns)
        result[arity] = np.array(weights)
    return result
