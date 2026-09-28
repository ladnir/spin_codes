"""Exact averaged support counts for independently permuted outer rows.

Separate ensemble from the parent directory's shared-column route.
No inner bound or full-code certificate is asserted here.
"""
from fractions import Fraction as Q
from itertools import accumulate
from math import comb, lcm, prod


def union_shells(spectra):
    """Expected tuple counts by union size, including the all-zero tuple.

    Each nonnegative integer spectrum assigns equal mass A[w]/C(n,w)
    to every support of size w. Independent row measures are multiplied.
    Exact integer differences avoid cancellation error in subset inversion.
    Componentwise spectrum uppers give componentwise *shell* uppers.
    """
    spectra = [tuple(row) for row in spectra]
    if not spectra or not spectra[0]:
        raise ValueError('at least one nonempty spectrum required')
    n = len(spectra[0]) - 1
    if any(len(row) != n+1 or any(not isinstance(x, int) or x < 0 for x in row)
           for row in spectra):
        raise ValueError('equal-length nonnegative integer spectra required')
    denominator = lcm(*(comb(n, w) for w in range(n+1)
                        if any(row[w] for row in spectra)))
    weighted = [[row[w] * (denominator // comb(n, w)) if row[w] else 0
                 for w in range(n+1)] for row in spectra]
    contained = [prod(sum(a * comb(v, w) for w, a in enumerate(row[:v+1]))
                      for row in weighted) for v in range(n+1)]
    divisor = denominator**len(spectra)
    result = []
    for u in range(n+1):
        numerator = comb(n, u) * contained[0]
        assert numerator >= 0, 'exact subset inversion must be nonnegative'
        result.append(Q(numerator, divisor))
        contained = [b-a for a, b in zip(contained, contained[1:])]
    assert sum(result) == prod(sum(row) for row in spectra)
    return result


def lowest_weight_spectrum(caps, total):
    """CDF-dominating spectrum with known total, NOT shellwise domination.

    Greedily spend the known total at the smallest permitted weights.
    Uniform subsets can be nested as weight increases. Consequently the
    union CDF from independent copies dominates the true union CDF.
    """
    if total < 0 or any(x < 0 for x in caps) or sum(caps) < total:
        raise ValueError('inconsistent spectrum cap or total')
    remaining = total
    result = []
    for cap in caps:
        value = min(cap, remaining)
        result.append(value)
        remaining -= value
    assert remaining == 0
    return result


def union_cdf_upper(caps, total, rows=4):
    """Exact rational CDF upper; excludes the all-zero tuple.

    Assumes A[0] = caps[0] = 1 and total = sum A. The returned CDF may be
    used by summation by parts, not by treating its differences as shell caps.
    """
    if not caps or caps[0] != 1 or rows < 1 or total < 1:
        raise ValueError('one zero word and a positive row count required')
    synthetic = lowest_weight_spectrum(caps, total)
    shells = union_shells([synthetic]*rows)
    shells[0] -= 1
    result = list(accumulate(shells))
    assert result[0] == 0 and result[-1] == total**rows-1
    return result


def fixed_weight_union(n, weights):
    if n < 0 or not weights or any(not 0 <= w <= n for w in weights):
        raise ValueError('invalid row weights')
    spectra = [[int(v == w) for v in range(n+1)] for w in weights]
    return union_shells(spectra)


def weighted_union_shells(spectrum, rows=4, full_weight=Q(1)):
    """Expected counts weighted by full_weight**J, J = all-row intersection.

    Full_weight >= 1 supports the old inner envelope's reciprocal all-one
    penalty. This uses only a supplied nonnegative spectrum; shell caps can
    replace the exact spectrum because the original subset sum is positive.
    """
    full_weight = Q(full_weight)
    if full_weight < 1 or rows < 1 or not spectrum:
        raise ValueError('positive row count and full_weight >= 1 required')
    if any(not isinstance(x, int) or x < 0 for x in spectrum):
        raise ValueError('nonnegative integer spectrum required')
    if full_weight == 1:
        return union_shells([spectrum]*rows)
    n = len(spectrum)-1
    denominator = lcm(*(comb(n,w) for w,a in enumerate(spectrum) if a))
    scaled = [a*(denominator//comb(n,w)) if a else 0 for w,a in enumerate(spectrum)]
    extra = full_weight-1
    a, b = extra.numerator, extra.denominator
    factors = [a**j * b**(n-j) for j in range(n+1)]
    contained = []
    for u in range(n+1):
        value = 0
        for j in range(u+1):
            row = sum(scaled[w]*comb(u-j,w-j) for w in range(j,u+1))
            value += comb(u,j)*factors[j]*row**rows
        contained.append(value)
    divisor = denominator**rows * b**n
    shells = []
    for u in range(n+1):
        assert contained[0] >= 0
        shells.append(Q(comb(n,u)*contained[0],divisor))
        contained = [b-a for a,b in zip(contained,contained[1:])]
    return shells


def weighted_cdf_upper(caps, total, rows=4, full_weight=Q(1)):
    """Nonzero-tuple weighted CDF; exact rationals, not shellwise caps."""
    full_weight = Q(full_weight)
    unweighted = union_cdf_upper(caps,total,rows)
    if full_weight == 1:
        return unweighted
    shells = weighted_union_shells(caps,rows,full_weight)
    shells[0] -= 1
    result = [min(v,unweighted[u]*full_weight**u)
              for u,v in enumerate(accumulate(shells))]
    assert all(a <= b for a,b in zip(result,result[1:]))
    return result


def expected_union(n, weights):
    if n <= 0:
        raise ValueError('positive length required')
    return n * (1-prod(Q(n-w, n) for w in weights))


def integer_cdf(values):
    """Upward integer interface for the existing outward CDF folder."""
    if any(x < 0 for x in values) or any(a>b for a,b in zip(values,values[1:])):
        raise ValueError('nonnegative nondecreasing CDF required')
    result=[-(-Q(x).numerator//Q(x).denominator) for x in values]
    assert all(a <= b for a,b in zip(values,result))
    assert all(a <= b for a,b in zip(result,result[1:]))
    return result


def tilted_union_shells(spectrum, rows=4, full_weight=Q(1), input_weight=Q(1)):
    """Exact shell counts weighted by input_weight**(-W) * full_weight**J.

    W is the total weight of the independently shuffled rows; J is the size
    of their common intersection. The all-zero tuple is included. Spectrum
    entries may be nonnegative rationals. Scaling each row measure to
    integer coefficients permits reuse of the exact subset inversion.

    Componentwise spectrum caps remain valid shell caps: the original
    expectation sums products of nonnegative row masses.
    """
    full_weight, input_weight = Q(full_weight), Q(input_weight)
    spectrum = tuple(Q(x) for x in spectrum)
    if (not spectrum or any(x < 0 for x in spectrum)
            or not isinstance(rows, int) or rows < 1
            or full_weight < 1 or input_weight <= 0):
        raise ValueError('nonnegative spectrum, positive integer rows/input weight, and full_weight >=1 required')
    tilted = [a / input_weight**w for w,a in enumerate(spectrum)]
    denominator = lcm(*(a.denominator for a in tilted))
    integral = [a.numerator * (denominator // a.denominator) for a in tilted]
    divisor = denominator**rows
    return [value / divisor for value in weighted_union_shells(integral,rows,full_weight)]


def tilted_support_caps(caps, total, rows=4, full_weight=Q(1), input_weight=Q(1)):
    """Return nonzero-tuple (shell caps, CDF caps) for the joint weight.

    The unknown spectrum A satisfies A[0]=caps[0]=1, A<=caps, and sum A=total.
    Shell caps use the componentwise upper spectrum. For input_weight >=1,
    the CDF also uses the known total: the lowest-weight spectrum dominates
    E[input_weight**(-W) * 1{union size <= u}] by a nested-subset coupling.
    Multiplying this bound by full_weight**u accounts for J<=u. For input
    weight below one, use W<=rows*u and J<=u with the unweighted known-total
    CDF instead; the decreasing-function coupling no longer applies.

    The product input_weight**(-W) * full_weight**J need not be monotone in
    the row weights. Consequently the greedy spectrum is not substituted
    directly into the joint weighted expectation.
    """
    caps = tuple(caps)
    full_weight, input_weight = Q(full_weight), Q(input_weight)
    if (not caps or caps[0] != 1
            or any(not isinstance(x, int) or x < 0 for x in caps)
            or not isinstance(total, int) or total < 1
            or not isinstance(rows, int) or rows < 1
            or full_weight < 1 or input_weight <= 0):
        raise ValueError('integer spectrum caps, one zero word, positive total/rows/input weight, and full_weight >=1 required')
    synthetic = lowest_weight_spectrum(caps,total)
    shells = tilted_union_shells(caps,rows,full_weight,input_weight)
    shells[0] -= 1
    if input_weight >= 1:
        decreasing_shells = tilted_union_shells(synthetic,rows,Q(1),input_weight)
        decreasing_shells[0] -= 1
        total_bounds = [decreasing * full_weight**u
                        for u,decreasing in enumerate(accumulate(decreasing_shells))]
    else:
        unweighted = union_cdf_upper(caps,total,rows)
        factor = full_weight / input_weight**rows
        total_bounds = [count * factor**u for u,count in enumerate(unweighted)]
    cdf = [min(joint,total_bound) for joint,total_bound in zip(accumulate(shells),total_bounds)]
    assert shells[0] == cdf[0] == 0
    assert all(x >= 0 for x in shells)
    assert all(a <= b for a,b in zip(cdf,cdf[1:]))
    return shells,cdf
