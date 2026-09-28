"""Prototype mass-based mature density bounds; not imported by verify.py.

For a fixed input shape, let beta=max_r Pr[BX=r], including r=0.
For each source expansion weight v, tilted feedback to any target is
at most beta*exp(-lambda*|v-W|). This gives a bound using mature mass
and its two nested tail masses. See ../MASS_DENSITY.md.

The alternative replaces the complete mature contribution to C, not
individual entries of the old column. Fixed convex mixtures are valid;
entrywise minima of these coupled alternatives generally are not.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb
from numbers import Integral, Rational
from pathlib import Path
import sys

from flint import arb

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from group_rank_one_verify import up
from occupancy_memory import Z, M, C
from mature_tail import L48, L56


LEVELS = (48, 56, 64, 72, 80)


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (Rational, str)):
        raise ValueError('exact rational parameters required')
    return Q(value)


def point(value):
    value = rational(value)
    return arb(value.numerator)/value.denominator


def upper(value):
    return up(value) if isinstance(value, arb) else value


def nested_coefficients(values, lower=48, middle=56):
    """Coefficients of M,L_middle,L_lower; rational or outward Arb inputs.

Each value must already upper-bound its expansion-class moment.
The three nonnegative coefficients dominate every supplied class.
"""
    if lower >= middle or not values or any(v < 0 for v in values.values()):
        raise ValueError('ordered cutoffs and nonnegative class moments required')
    values = {v:upper(value) for v,value in values.items()}
    a = max((f for v,f in values.items() if v > middle), default=0)
    b = max(0, upper(max((f for v,f in values.items() if lower < v <= middle), default=0)-a))
    c = max(0, upper(max((f for v,f in values.items() if v <= lower), default=0)-a-b))
    return a,b,c


def coefficients(data, maximum, tilt, penalty='1', input_penalty='1', *, rounds=2):
    """Complete shape maxima using full_feedback_census.census() records.

The caller must generate these records for the same authenticated maps.
This function checks completeness and elementary exact counting guards;
it does not authenticate an arbitrary external census file.
"""
    if not isinstance(maximum, int) or isinstance(maximum, bool) or not 1 <= maximum <= 32:
        raise ValueError('one through 32 packet windows required')
    if not isinstance(rounds, int) or isinstance(rounds, bool) or rounds < 1:
        raise ValueError('positive integer update count required')
    tilt, penalty, input_penalty = map(rational,(tilt,penalty,input_penalty))
    if tilt < 0 or not 0 < penalty <= 1 or input_penalty <= 0:
        raise ValueError('require lambda>=0, 0<rho<=1, and input weight factor>0')
    required = {shape for j in range(1,maximum+1)
                for shape in combinations_with_replacement(range(1,5),j)}
    if not required <= set(data):
        raise ValueError('complete feedback census required through requested occupancy')
    peaks = {}
    for shape in sorted(required):
        zero,peak,denominator,_ = data[shape]
        if (any(not isinstance(v, Integral) or isinstance(v, bool) for v in (zero,peak,denominator))
                or denominator <= 0 or min(zero,peak) < 0 or max(zero,peak) > denominator):
            raise ValueError('invalid exact feedback counts')
        peaks[shape] = Q(max(zero,peak),denominator)
    return _coefficients_from_peaks(peaks,maximum,tilt,penalty,input_penalty,rounds)


def _coefficients_from_peaks(peaks, maximum, tilt, penalty, input_penalty, rounds):
    """Internal complete-shape calculation after the caller's validation."""
    powers = {w:upper((-point(tilt)*w).exp()) for w in range(129)}
    result = {j:[arb(0),arb(0),arb(0)] for j in range(1,maximum+1)}
    for shape,beta in peaks.items():
        scale = point(beta*penalty**shape.count(4)*input_penalty**sum(shape)/2**rounds)
        local = nested_coefficients({v:powers[abs(v-sum(shape))] for v in LEVELS})
        for i,value in enumerate(local):
            result[len(shape)][i] = max(result[len(shape)][i], upper(scale*value))
    return {j:tuple(row) for j,row in result.items()}


def conditioned_peaks(data, through, maximum, *, windows=32, packet_bits=4):
    """Complete shape bounds using a smaller exact feedback census.

Expose j-m labeled packets first. The remaining m packets occupy uniform
distinct windows among windows-j+m positions. Conditioning an unrestricted
m-packet injection to avoid the exposed positions costs at most
binom(windows,m)/binom(windows-j+m,m). Shifting its feedback does not
change its maximum atom. Minimize over all weight submultisets available
in the j-packet shape, always including the zero atom in census peaks.

The dynamic program finds these submultiset minima without enumerating
every pair of full shapes and subshapes. All values are exact rationals.
"""
    if (any(not isinstance(v,int) or isinstance(v,bool) for v in (through,maximum,windows,packet_bits))
            or not 1 <= through <= maximum <= windows or packet_bits < 1):
        raise ValueError('require 1 <= census cutoff <= maximum <= windows and positive packet width')
    exact = {}
    for j in range(1,through+1):
        for shape in combinations_with_replacement(range(1,packet_bits+1),j):
            if shape not in data:
                raise ValueError('complete feedback census required through cutoff')
            zero,peak,denominator,_ = data[shape]
            if (any(not isinstance(v,Integral) or isinstance(v,bool) for v in (zero,peak,denominator))
                    or denominator <= 0 or min(zero,peak)<0 or max(zero,peak)>denominator):
                raise ValueError('invalid exact feedback counts')
            exact[shape] = Q(max(zero,peak),denominator)
    result, previous = {},{}
    for j in range(1,maximum+1):
        current = {m:{} for m in range(1,min(j,through)+1)}
        factors = {m:Q(comb(windows,m),comb(windows-j+m,m)) for m in current}
        for shape in combinations_with_replacement(range(1,packet_bits+1),j):
            children = [shape[:i]+shape[i+1:] for i,b in enumerate(shape) if i==0 or b!=shape[i-1]]
            for m in current:
                current[m][shape] = exact[shape] if j==m else min(previous[m][child] for child in children)
            result[shape] = exact[shape] if j<=through else min(Q(1),*(current[m][shape]*factors[m] for m in current))
        previous = current
    return result


def conditioned_coefficients(data, through, maximum, tilt, penalty='1', input_penalty='1', *, rounds=2):
    """Mass-density columns through maximum using a census through cutoff."""
    # Reuse validation of the outward parameters and exact census records.
    coefficients(data,through,tilt,penalty,input_penalty,rounds=rounds)
    peaks = conditioned_peaks(data,through,maximum)
    tilt,penalty,input_penalty = map(rational,(tilt,penalty,input_penalty))
    return _coefficients_from_peaks(peaks,maximum,tilt,penalty,input_penalty,rounds)


def blend(base, candidates, fractions, *, target=C):
    """Fixed per-occupancy mass fractions in [0,1]; keep other columns.

Fractions must be fixed proof witnesses, not functions of incoming state.
Apply only after the existing unsplit density refinements. No default
heuristic or proof claim chooses these fractions for the caller.

For target=Z, the baseline C->Z entry must bound the lazy mature return;
its M/L->Z entries must independently bound the refresh contribution.
Retain those refresh coefficients and add the alternative lazy column.
The same peak bound includes the zero target. This option is specific to
the baseline's stated component decomposition, not arbitrary matrices.
"""
    if target not in (Z,C):
        raise ValueError('only the mature-density and zero-return columns are supported')
    if any(t.nrows() != 11 or t.ncols() != 11 for t in base):
        raise ValueError('eleven-coordinate operators required')
    result = [t*1 for t in base]
    for j,fraction in fractions.items():
        if (not isinstance(j,int) or isinstance(j,bool) or not 1 <= j < len(base)
                or j not in candidates):
            raise ValueError('nonzero covered local occupancy required')
        fraction = rational(fraction)
        if not 0 <= fraction <= 1:
            raise ValueError('convex fraction outside [0,1]')
        if any(base[j][source,C] != 0 for source in (M,L56,L48)):
            raise ValueError('input mature-density contribution must be unsplit')
        row = candidates[j]
        if len(row) != 3 or any(not x.is_finite() or not x >= 0 for x in row):
            raise ValueError('three finite nonnegative outward coefficients required')
        if fraction == 0:
            continue
        for source,value in zip((M,L56,L48),row):
            refresh = base[j][source,Z] if target==Z else 0
            result[j][source,target] = upper(refresh+point(fraction)*value)
        result[j][C,target] = upper(point(1-fraction)*base[j][C,target])
    return result
