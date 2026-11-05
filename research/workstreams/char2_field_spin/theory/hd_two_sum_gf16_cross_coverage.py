"""Positive scalar contraction requiring every cross-state region active.

No regional provider, fitting, production driver, or certificate. A_j is
the weight of ONE FIXED j-subset of c labelled cross regions after all
optional coordinates have been summed. It does not include C(c,j).
The group law must be symmetric in those c coordinates. Group/type laws
are identical only in the comparison; no actual row sets are resampled.
"""
from functools import lru_cache
from math import comb, log2

import numpy as np


@lru_cache(None)
def coverage_geometry(cross):
    """Immutable (oldcovered,newcovered,[(size,log multiplicity)]) rows."""
    if type(cross) is not int or not 0 <= cross <= 16:
        raise ValueError('The declared cross-region count must be in0..16')
    return tuple((old, old+new,
                  tuple((old_used+new, log2(comb(old, old_used)*comb(cross-old, new)))
                        for old_used in range(old+1)))
                 for old in range(cross+1) for new in range(cross-old+1))


def subset_weights(weights, active_zero, cross, live, log_zero_atom,
                   log_cross_atom, *, log_marker=0., arity=8):
    """Return log A_j for a FIXED cross j-subset, with scalar floor mark.

    A_j=Bcross^j sum_(u,v) C(k,u)C(l,v)W_(u+j+v)Bzero^u
                                  marker^min(u+j+v,arity).
    Live-state atom is1. Remaining coordinates are forced empty.
    No cross-subset binomial belongs in this function.
    """
    values = np.asarray(weights, dtype=float)
    if (values.ndim != 1 or len(values) < 2 or len(values) > 17
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))
            or any(type(x) is not int or x < 0 for x in (active_zero, cross, live))
            or active_zero+cross+live >= len(values)
            or type(arity) is not int or not 1 <= arity < len(values)
            or any(not np.isfinite(x) for x in (log_zero_atom, log_cross_atom, log_marker))):
        raise ValueError('Valid labelled coordinate counts, weights, and finite scalar marks required')
    coverage_geometry(cross)
    answer = np.full(cross+1, -np.inf)
    for selected in range(cross+1):
        terms = [log2(comb(active_zero, u)*comb(live, v))
                 +values[u+selected+v]+u*log_zero_atom+selected*log_cross_atom
                 +min(u+selected+v, arity)*log_marker
                 for u in range(active_zero+1) for v in range(live+1)]
        answer[selected] = np.logaddexp2.reduce(terms)
    return answer


def coverage_transition(log_subset_weights):
    """Positive covered-count transition; no inclusion-exclusion."""
    values = np.asarray(log_subset_weights, dtype=float)
    if (values.ndim != 1 or not 1 <= len(values) <= 17
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))):
        raise ValueError('Complete A0..Ac finite-or-zero fixed-subset weights required')
    cross = len(values)-1
    matrix = np.full((cross+1, cross+1), -np.inf)
    for old, new, terms in coverage_geometry(cross):
        matrix[old, new] = np.logaddexp2.reduce([values[size]+coefficient for size, coefficient in terms])
    return matrix


def covered_groups(log_subset_weights, groups):
    """log(sum product A_|Sg|) over unions S1∪...∪Sgroups=[c]."""
    if type(groups) is not int or groups < 0:
        raise ValueError('A nonnegative integer number of groups is required')
    transition = coverage_transition(log_subset_weights)
    cross = len(transition)-1
    values = np.full(cross+1, -np.inf)
    values[0] = 0.
    for _ in range(groups):
        values = np.logaddexp2.reduce(values[:, None]+transition, axis=0)
    return float(values[cross])


def unrestricted_groups(log_subset_weights, groups):
    """Matched no-coverage control, including C(c,j) once."""
    coverage_transition(log_subset_weights)
    if type(groups) is not int or groups < 0:
        raise ValueError('A nonnegative integer number of groups is required')
    if groups == 0:
        return 0.
    cross = len(log_subset_weights)-1
    return groups*float(np.logaddexp2.reduce([
        value+log2(comb(cross, selected))
        for selected, value in enumerate(log_subset_weights)]))
