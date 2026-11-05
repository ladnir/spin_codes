"""Bound a common-potential support sum by cyclic runs and convolution.

Fix scalar support weight r and one positive state vector v, common to all
matching histograms. Write c_f for a contraction of the regional matrix
with f full pairs, and p_f=c_f**128. Every indexed support with matching
counts a,b has r-a-b cyclic runs. Hence its class count is at most the
whole-cycle count C(2M,r,r-a-b). This gives

  sum_U p_a(U) p_b(U) <= sum_d C(2M,r,r-d) [X**d] P(X)**2,
  P(X) = sum_f p_f X**f.

The coefficient product does not assert independent matchings. It sums
upper bounds on every joint class. In particular the impossible zero-run
class vanishes for nonconstant supports. Multiplying by v[0]/min(v)
bounds the original ordered regional matrix product. This is a bound on
one weight class, not a whole-code distance certificate.
"""
from functools import lru_cache
from math import log2

import numpy as np
from flint import arb,arb_poly

from hd_shifted_alias_certificate import rational
from hd_shifted_pairs_count import circular_run_count,one_pairing_count
from hd_shifted_swap_potential import full_counts,positive_vector,contraction_arb


def circular_run_counts(length, weight):
    """Exact indexed-cycle counts by run count, using integer recurrence.

    Nonconstant cycles have N[1]=length and
    N[c+1]=N[c]*(weight-c)*(length-weight-c)/(c*(c+1)).
    This is the ratio of the same binomial formula used by
    circular_run_count. Every division is checked for exactness.
    """
    if (type(length) is not int or type(weight) is not int
            or length < 1 or not 0 <= weight <= length):
        raise ValueError('An integer cycle length and admissible weight are required')
    counts = [0]*(weight+1)
    if weight in (0, length):
        counts[0] = 1
        return counts
    counts[1] = length
    for runs in range(1, min(weight, length-weight)):
        numerator = counts[runs]*(weight-runs)*(length-weight-runs)
        value, remainder = divmod(numerator, runs*(runs+1))
        if remainder:
            raise ArithmeticError('Nonintegral indexed-cycle count recurrence')
        counts[runs+1] = value
    return counts


def run_sum_arb(pairs,weight,regions,vector,*,repetitions=256):
    """Directed positive polynomial convolution at the caller's precision."""
    vector = positive_vector(vector)
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('Positive even regional repetition count required')
    counts = full_counts(pairs,weight)
    coefficients = [arb(0) for _ in range(weight//2+1)]
    for full in counts:
        coefficients[full] = (contraction_arb(regions[full],vector)**(repetitions//2)).upper()
    convolution = arb_poly(coefficients)**2
    counts_by_runs = circular_run_counts(2*pairs,weight)
    value = arb(0)
    for degree in range(2*counts.start,2*(counts.stop-1)+1):
        multiplicity = counts_by_runs[weight-degree]
        if not multiplicity:
            continue
        upper = convolution[degree].upper()
        if not upper>=0 or not upper.is_finite():
            raise ArithmeticError('Invalid directed positive convolution coefficient')
        value += multiplicity*upper
    return (value*rational(vector[0]/min(vector))).upper()


@lru_cache(8)
def log_run_counts(pairs,weight):
    counts = tuple(full_counts(pairs,weight))
    degrees = np.add.outer(counts,counts)
    row = []
    for degree in range(weight+1):
        value = circular_run_count(2*pairs,weight,weight-degree)
        row.append(log2(value) if value else -np.inf)
    result = np.asarray(row)[degrees]
    result.setflags(write=False)
    return result


def run_log_screen(pairs,weight,log_regions,log_ratio,*,repetitions=256):
    """FLOAT proposal only; direct log sums avoid unstable FFT subtraction."""
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('Positive even regional repetition count required')
    entries = np.asarray([np.asarray(log_regions[f]).ravel() for f in full_counts(pairs,weight)])
    a,b,c,d = entries.T
    logs = repetitions/2*np.maximum(np.logaddexp2(a,b+log_ratio),
                                    np.logaddexp2(c-log_ratio,d))
    terms = log_run_counts(pairs,weight)+np.add.outer(logs,logs)
    return float(np.logaddexp2.reduce(terms.ravel())+max(0.,-log_ratio))


@lru_cache(4)
def log_joint_counts(pairs,weight):
    """FLOAT version of the established profile_count_upper, for diagnosis."""
    counts = tuple(full_counts(pairs,weight))
    marginal = np.array([log2(one_pairing_count(pairs,weight,f)) for f in counts])
    result = np.minimum(log_run_counts(pairs,weight),np.minimum.outer(marginal,marginal))
    if weight%2==0:
        # One all-full matching contracts the cycle from2M toM positions.
        half = weight//2
        endpoint = []
        for f in counts:
            value = circular_run_count(pairs,half,half-f)
            endpoint.append(log2(value) if value else -np.inf)
        result[-1,:] = endpoint
        result[:,-1] = endpoint
    result.setflags(write=False)
    return result


def joint_log_screen(pairs,weight,log_regions,log_ratio,*,repetitions=256):
    """FLOAT diagnostic separating run-count slack from potential loss."""
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('Positive even regional repetition count required')
    entries = np.asarray([np.asarray(log_regions[f]).ravel() for f in full_counts(pairs,weight)])
    a,b,c,d = entries.T
    logs = repetitions/2*np.maximum(np.logaddexp2(a,b+log_ratio),
                                    np.logaddexp2(c-log_ratio,d))
    terms = log_joint_counts(pairs,weight)+np.add.outer(logs,logs)
    return float(np.logaddexp2.reduce(terms.ravel())+max(0.,-log_ratio))
