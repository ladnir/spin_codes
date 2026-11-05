"""Directed scalar envelopes for boxes of alternating-route support profiles.

This is one component of a future coverage certificate, not a code bound.
Regional coefficient matrices and support counts must be supplied separately.
Convex log-factorial chords retain the shared profile dependence of both
regional normalizers instead of maximizing their reciprocals separately.
"""

from fractions import Fraction
from itertools import product
from math import factorial

from flint import arb, ctx

from hd_shifted_general_count import ProfileBox, matching_histograms


def histogram_affines(pairs):
    """Each entry is (constant, coefficients of E,O,a,b)."""
    return (
        ((pairs,(-1,-1,1,0)), (0,(1,0,-1,0)),
         (0,(0,1,-1,0)), (0,(0,0,1,0))),
        ((pairs,(-1,-1,0,1)), (0,(0,1,0,-1)),
         (0,(1,0,0,-1)), (0,(0,0,0,1))),
    )


def rational(value):
    value = Fraction(value)
    return arb(value.numerator)/value.denominator


def histogram_ranges(box,pairs=4096):
    """Integer envelopes of both oriented histograms over a profile box."""
    result = []
    for histogram in histogram_affines(pairs):
        ranges = []
        for offset,coefficients in histogram:
            lo = offset+sum(c*(lo if c>=0 else hi) for c,(lo,hi) in zip(coefficients,box.ranges))
            hi = offset+sum(c*(hi if c>=0 else lo) for c,(lo,hi) in zip(coefficients,box.ranges))
            lo,hi = max(0,lo),min(pairs,hi)
            if lo>hi:
                raise ValueError('Empty histogram envelope')
            ranges.append((lo,hi))
        result.append(tuple(ranges))
    return tuple(result)


def scalar_envelope(pairs, box, z, alpha, markers, *, repetitions=128,
                    message_symbols_per_group=128, weight_limit=734003,
                    q=256, precision=256):
    """Upper-bound the profile-dependent logarithmic scalar over a box.

    The bounded quantity is
      alpha*(r*(2*repetitions*log2(q)-k*log2(q-1))
             -log2(q-1)-weight_limit*log2(z))
      -repetitions*sum_j log2(multinomial(H_j)*markers_j**H_j[1:]).

    A regional alias matrix is unnormalized here; its marker monomial and
    multinomial have both been placed in this scalar. No operator envelope
    or class multiplicity is included. The bound applies to all feasible
    profiles in the box even when some of its corners are infeasible.
    """
    if not isinstance(pairs,int) or pairs<1 or not isinstance(box,ProfileBox):
        raise ValueError('Expected a positive pair count and ProfileBox')
    if any(len(interval)!=2 or any(not isinstance(n,int) or not 0<=n<=pairs for n in interval)
           or interval[0]>interval[1] for interval in box.ranges):
        raise ValueError('Invalid profile box')
    z,alpha = Fraction(z),Fraction(alpha)
    markers = tuple(tuple(map(Fraction,row)) for row in markers)
    if (not 0<z<1 or not 0<alpha<=1 or len(markers)!=2
            or any(len(row)!=3 or min(row)<=0 for row in markers)
            or not isinstance(repetitions,int) or repetitions<1
            or not isinstance(message_symbols_per_group,int) or message_symbols_per_group<1
            or not isinstance(weight_limit,int) or weight_limit<0
            or not isinstance(q,int) or q<2 or precision<64):
        raise ValueError('Invalid scalar-envelope parameters')
    with ctx.workprec(precision):
        log2 = lambda value:value.log()/arb(2).log()
        logs = {0:arb(0),1:arb(0)}
        def log_factorial(n):
            if n not in logs:
                logs[n] = log2(arb(factorial(n)))
            return logs[n]
        density = 2*repetitions*log2(arb(q))-message_symbols_per_group*log2(arb(q-1))
        exponent = rational(alpha)
        constant = -exponent*(log2(arb(q-1))+weight_limit*log2(rational(z)))
        constant -= 2*repetitions*log_factorial(pairs)
        slopes = [exponent*density,exponent*density,arb(0),arb(0)]
        for parity,histogram in enumerate(histogram_affines(pairs)):
            for index,(offset,coefficients) in enumerate(histogram):
                lower = offset+sum(c*(lo if c>=0 else hi)
                                   for c,(lo,hi) in zip(coefficients,box.ranges))
                upper = offset+sum(c*(hi if c>=0 else lo)
                                   for c,(lo,hi) in zip(coefficients,box.ranges))
                # Only feasible profiles matter; their histogram counts
                # lie in [0,pairs], even when the box also contains others.
                lower,upper = max(0,lower),min(pairs,upper)
                if lower>upper:
                    raise ValueError('Box contains no feasible histogram')
                if lower==upper:
                    chord_slope = arb(0)
                    chord_intercept = log_factorial(lower)
                else:
                    chord_slope = (log_factorial(upper)-log_factorial(lower))/(upper-lower)
                    chord_intercept = log_factorial(lower)-lower*chord_slope
                term_slope = repetitions*chord_slope
                if index:
                    term_slope -= repetitions*log2(rational(markers[parity][index-1]))
                constant += repetitions*chord_intercept+offset*term_slope
                for coordinate,c in enumerate(coefficients):
                    slopes[coordinate] += c*term_slope
        # Each affine candidate is enclosed. Selecting the largest upper
        # endpoint is directed even when an interval slope straddles zero.
        corners = product(*(tuple(dict.fromkeys(interval)) for interval in box.ranges))
        upper = max((constant+sum((s*x for s,x in zip(slopes,corner)),arb(0))).upper()
                    for corner in corners)
        return upper


def point_scalar(pairs, profile, z, alpha, markers, *, repetitions=128,
                 message_symbols_per_group=128, weight_limit=734003,
                 q=256, precision=256):
    """Independent direct formula for testing the chord envelope."""
    histograms = matching_histograms(pairs,*profile)
    with ctx.workprec(precision):
        log2 = lambda value:value.log()/arb(2).log()
        value = rational(alpha)*(sum(profile[:2])*(2*repetitions*log2(arb(q))
            -message_symbols_per_group*log2(arb(q-1)))-log2(arb(q-1))-weight_limit*log2(rational(z)))
        for histogram,row in zip(histograms,markers):
            normalizer = Fraction(factorial(pairs))
            for count in histogram:
                normalizer /= factorial(count)
            for count,marker in zip(histogram[1:],row):
                normalizer *= Fraction(marker)**count
            value -= repetitions*log2(rational(normalizer))
        return value
