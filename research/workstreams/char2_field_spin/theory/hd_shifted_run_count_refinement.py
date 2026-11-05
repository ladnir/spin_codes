"""Count refinements for a fixed scalar support weight; no encoder changes.

For d=a+b, the joint classes partition C(2M,r,r-d) indexed supports.
Thus sum_(a+b=d) N_ab p_a p_b <= C_d max_(a+b=d) p_a p_b.
For any positive integer q, replace the maximum by
([X^d](sum_f p_f^q X^f)^2)^(1/q). This needs one positive convolution,
not an independence assumption. q=1 is the existing run bound.

An exact joint-count diagnostic uses alternating Markov types. Fix E ones
on even positions and O=r-E on odd positions. The even-edge table is
(M-r+a,O-a,E-a,a); the odd-edge table is (M-r+b,E-b,O-b,b).
The directed multigraph has vertices (even0,even1,odd0,odd1), outdegrees
(M-E,E,M-O,O), and these two tables as its edge multiplicities. Even
edges run from even to odd vertices; odd edges run back to even vertices.
For four occupied bit/parity states, a last-exit-tree count gives

 N_E = M*tau * product_v (degree_v-1)! / product_edges count_edge!,
 tau = E*O*(M-r+a+b)-M*a*b.

The outgoing Laplacian minor with root even0 is

 [[E, -(E-a), -a], [-(E-b), M-O, 0], [-b, 0, O]].

Its determinant is tau. The balanced graph has the same in-tree cofactor
for each root; a disconnected type has tau=0 and contributes no word.
To justify the factor: count closed edge sequences rooted at an even
state. Last exits from the other states form an in-tree, and the remaining
outgoing edges can be ordered arbitrarily. Divide by permutations of
identical transitions. The root contributes degree_root!, and each other
vertex contributes (degree-1)!, so summing roots even0/even1 contributes
(M-E)+E=M. Position zero has fixed even parity; no further rotation or
starting-position quotient is taken. If E=0 or O=0, a=b=0 and the count is
binom(M,max(E,O)). If E=M or O=M, a=b=min(E,O) and the count is
binom(M,min(E,O)). These rules also cover the constant words. Summing N_E
over E gives the indexed joint class count. Exact toy tests independently
enumerate all oriented (r,a,b,E) classes.

The exact count is a diagnostic, not yet a fast replacement for every
profile. The convolution refinement reuses a fixed common potential and
does not itself include the outer/projective prefix.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, lgamma, log, log2

import numpy as np
from flint import arb, arb_poly

from hd_shifted_alias_certificate import rational
from hd_shifted_pairs_count import circular_run_count
from hd_shifted_swap_potential import contraction_arb, full_counts, positive_vector


@lru_cache(2)
def _factorials(pairs):
    values = [1]
    for n in range(1, pairs+1):
        values.append(values[-1]*n)
    return tuple(values)


def conditional_profile_count(pairs, weight, a, b, even):
    """Exact number of indexed supports with the specified even weight."""
    if (any(not isinstance(x, int) for x in (pairs, weight, a, b, even))
            or pairs < 1):
        raise ValueError('Integer parameters and a positive pair count required')
    odd = weight-even
    if not 0<=even<=pairs or not 0<=odd<=pairs:
        return 0
    edges = (pairs-weight+a, odd-a, even-a, a,
             pairs-weight+b, even-b, odd-b, b)
    if min(edges)<0:
        return 0
    if even==0 or odd==0:
        return comb(pairs, max(even, odd)) if a==b==0 else 0
    if even==pairs or odd==pairs:
        other = min(even, odd)
        return comb(pairs, other) if a==b==other else 0
    tau = even*odd*(pairs-weight+a+b)-pairs*a*b
    if tau==0:
        return 0
    if tau<0:
        raise ArithmeticError('Negative tree determinant for a valid type')
    facts = _factorials(pairs)
    numerator = pairs*tau
    for degree in (pairs-even, even, pairs-odd, odd):
        numerator *= facts[degree-1]
    denominator = 1
    for count in edges:
        denominator *= facts[count]
    answer, remainder = divmod(numerator, denominator)
    if remainder:
        raise ArithmeticError('Nonintegral indexed Markov-type count')
    return answer


def exact_profile_count(pairs, weight, a, b):
    """Exact joint class count via a one-dimensional even-weight sum."""
    if any(not isinstance(x, int) for x in (pairs, weight, a, b)) or pairs<1:
        raise ValueError('Integer parameters and a positive pair count required')
    if not 0<=weight<=2*pairs or min(a,b)<0 or max(a,b)>pairs:
        return 0
    low = max(0, weight-pairs, a, b)
    high = min(pairs, weight, weight-a, weight-b)
    return sum(conditional_profile_count(pairs,weight,a,b,even)
               for even in range(low,high+1))


def exact_joint_sum_arb(pairs,weight,regions,vector,*,repetitions=256):
    """Exact correlated count kernel with directed common-potential caps.

    For interior E,O, set h_a=binom(E,a)binom(M-E,O-a)/binom(M,O),
    V=binom(M,E)binom(M,O), S_j=sum_a h_a a^j p_a. The exact conditional
    sum is VM/((M-E)(M-O))*((M-r)S0^2+2S0*S1-M*S1^2/(EO)).
    This identity follows by factoring the last-exit count above. h_a*h_b
    is only an algebraic factor; the bracket retains their dependence.
    Arb encloses the subtraction, and only the final upper endpoint is
    used. Coefficientwise replacement of p by upper bounds is monotone
    because the underlying exact joint counts are nonnegative.
    """
    vector = positive_vector(vector)
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('Positive even repetitions required')
    counts = full_counts(pairs,weight)
    caps = {f:(contraction_arb(regions[f],vector)**(repetitions//2)).upper() for f in counts}
    value = arb(0)
    for even in range(max(0,weight-pairs),min(pairs,weight)+1):
        odd = weight-even
        if even==0 or odd==0:
            value += comb(pairs,max(even,odd))*caps[0]**2
            continue
        if even==pairs or odd==pairs:
            other = min(even,odd)
            value += comb(pairs,other)*caps[other]**2
            continue
        low,high = max(0,weight-pairs),min(even,odd)
        h = comb(even,low)*comb(pairs-even,odd-low)
        t0,t1 = arb(0),arb(0)
        for full in range(low,high+1):
            term = h*caps[full]
            t0 += term
            t1 += full*term
            if full<high:
                h,remainder = divmod(h*(even-full)*(odd-full),
                                     (full+1)*(pairs-weight+full+1))
                if remainder:
                    raise ArithmeticError('Nonintegral hypergeometric recurrence')
        bracket = (pairs-weight)*t0*t0+2*t0*t1-rational(F(pairs,even*odd))*t1*t1
        scale = rational(F(comb(pairs,even)*pairs,
                           comb(pairs,odd)*(pairs-even)*(pairs-odd)))
        upper = (scale*bracket).upper()
        if not upper>=0 or not upper.is_finite():
            raise ArithmeticError('Invalid exact-count moment enclosure')
        value += upper
    return (value*rational(vector[0]/min(vector))).upper()


@lru_cache(2)
def _log_factorials(pairs):
    return np.array([lgamma(n+1)/log(2) for n in range(pairs+1)])


def exact_joint_log_screen(pairs,weight,log_regions,log_ratio,*,repetitions=256):
    """FLOAT evaluation of the exact count kernel, using two moments per E."""
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('Positive even repetitions required')
    counts = full_counts(pairs,weight)
    entries = np.array([np.asarray(log_regions[f]).ravel() for f in counts])
    a,b,c,d = entries.T
    caps = repetitions/2*np.maximum(np.logaddexp2(a,b+log_ratio),np.logaddexp2(c-log_ratio,d))
    lf = _log_factorials(pairs)
    choose = lambda n,k:lf[n]-lf[k]-lf[n-k]
    terms = []
    for even in range(max(0,weight-pairs),min(pairs,weight)+1):
        odd = weight-even
        if even==0 or odd==0:
            terms.append(choose(pairs,max(even,odd))+2*caps[-counts.start])
            continue
        if even==pairs or odd==pairs:
            other = min(even,odd)
            terms.append(choose(pairs,other)+2*caps[other-counts.start])
            continue
        low,high = max(0,weight-pairs),min(even,odd)
        full = np.arange(low,high+1)
        logs = choose(even,full)+choose(pairs-even,odd-full)-choose(pairs,odd)+caps[full-counts.start]
        peak = np.max(logs)
        if not np.isfinite(peak):
            continue
        scaled = np.exp2(logs-peak)
        s0 = float(np.sum(scaled))
        mean = float(np.dot(full,scaled)/s0)
        mean = min(float(high),max(float(low),mean))
        # Concave-quadratic Bernstein form avoids subtracting near an endpoint.
        endpoint = lambda f:float(F(pairs-weight+2*f)-F(pairs*f*f,even*odd))
        fraction = (mean-low)/(high-low)
        bracket = ((1-fraction)*endpoint(low)+fraction*endpoint(high)
                   +pairs/(even*odd)*(mean-low)*(high-mean))
        if bracket<=0:
            continue
        terms.append(choose(pairs,even)+choose(pairs,odd)+log2(pairs)
                     -log2(pairs-even)-log2(pairs-odd)+2*(peak+log2(s0))+log2(bracket))
    return float(np.logaddexp2.reduce(terms)+max(0.,-log_ratio))


def run_partition_sum_arb(pairs, weight, regions, vector, *, q=2, repetitions=256):
    """Directed run-partition l_q bound at the caller's precision.

    One common positive vector applies to every histogram. q=1 recovers
    the existing convolution argument; q>1 reduces its count overcharge.
    """
    vector = positive_vector(vector)
    if (not isinstance(q,int) or q<1 or not isinstance(repetitions,int)
            or repetitions<2 or repetitions%2):
        raise ValueError('Positive integer q and positive even repetitions required')
    counts = full_counts(pairs,weight)
    coefficients = [arb(0) for _ in range(weight//2+1)]
    for full in counts:
        coefficients[full] = (contraction_arb(regions[full],vector)**(q*repetitions//2)).upper()
    convolution = arb_poly(coefficients)**2
    value = arb(0)
    for degree in range(2*counts.start,2*(counts.stop-1)+1):
        count = circular_run_count(2*pairs,weight,weight-degree)
        if not count:
            continue
        upper = convolution[degree].upper()
        if not upper>=0 or not upper.is_finite():
            raise ArithmeticError('Invalid positive convolution coefficient')
        term = upper if q==1 else (upper**rational(F(1,q))).upper()
        value += count*term
    return (value*rational(vector[0]/min(vector))).upper()


def run_partition_log_screen(pairs,weight,log_regions,log_ratio,*,q=2,repetitions=256):
    """Stable FLOAT diagnostic; q=None computes the exact maximum per run."""
    if q is not None and (not isinstance(q,int) or q<1):
        raise ValueError('Positive integer q or None required')
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('Positive even repetitions required')
    counts = full_counts(pairs,weight)
    entries = np.array([np.asarray(log_regions[f]).ravel() for f in counts])
    a,b,c,d = entries.T
    logs = repetitions/2*np.maximum(np.logaddexp2(a,b+log_ratio),np.logaddexp2(c-log_ratio,d))
    terms = []
    for degree in range(2*counts.start,2*(counts.stop-1)+1):
        count = circular_run_count(2*pairs,weight,weight-degree)
        if not count:
            continue
        lo,hi = max(counts.start,degree-counts.stop+1),min(counts.stop-1,degree-counts.start)
        full = np.arange(lo,hi+1)
        products = logs[full-counts.start]+logs[degree-full-counts.start]
        value = np.max(products) if q is None else np.logaddexp2.reduce(q*products)/q
        terms.append(log2(count)+value)
    return float(np.logaddexp2.reduce(terms)+max(0.,-log_ratio))
