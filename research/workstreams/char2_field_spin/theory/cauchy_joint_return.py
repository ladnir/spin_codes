"""Exact fixed-word joint return/output laws for one to three packet inputs.

The incoming state and nonzero input values are fixed BEFORE all current
packet diagonals.  This law consumes those diagonals jointly and must not
be multiplied into another law that has already conditioned on them.
No native implementation, benchmark, or whole-code certificate is supplied.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from math import comb

from cauchy_diagonal_packet_reference import LABELS, NODES
from rate_design_coset import mul, power


Q = 256
MUL = tuple(tuple(mul(a, b) for b in range(Q)) for a in range(Q))
INV = (0,)+tuple(power(a, Q-2) for a in range(1, Q))


def dot(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= MUL[x][y]
    return result


def canonical(vector):
    first = next((value for value in vector if value), 0)
    if not first:
        raise ValueError("zero vector has no projective representative")
    inverse = INV[first]
    return tuple(MUL[inverse][value] for value in vector)


def matrix_rank(matrix, columns):
    work = [list(row) for row in matrix]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = INV[work[rank][column]]
        work[rank] = [MUL[inverse][value] for value in work[rank]]
        for i in range(rank+1, len(work)):
            factor = work[i][column]
            work[i] = [x ^ MUL[factor][y] for x, y in zip(work[i], work[rank])]
        rank += 1
    return rank


def packet_forms(phase, inputs):
    if phase not in range(4) or not 1 <= len(inputs) <= 3 or tuple(sorted(set(inputs))) != tuple(inputs):
        raise ValueError("require phase0..3 and one to three distinct ordered input indices")
    if any(i not in range(16) for i in inputs):
        raise ValueError("input index outside packet")
    labels = LABELS[phase]
    injection = tuple(tuple(power(labels[i] ^ v, 3) for i in inputs) for v in NODES)
    returned_output = tuple(tuple(int(j == i) ^ MUL[2][power(labels[j] ^ labels[i], 3)]
                                  for i in inputs) for j in range(16))
    return injection, returned_output


def projective_strata(injection, returned_output):
    """Count projective inputs by (H, output-zero mask, full-zero rank).

    For dimension three this enumerates a projective-line arrangement:
    intersection points, points on exactly one line, and points on none.
    Its complexity depends on the number of forms, not on q^2 vectors.
    """
    w = len(injection[0])
    if w == 1:
        h = sum((row[0] != 0) << j for j, row in enumerate(injection))
        zeros = sum((row[0] == 0) << j for j, row in enumerate(returned_output))
        return {(h, zeros, 0): 1}
    if w not in (2, 3):
        raise ValueError("only projective dimensions one and two are supported")
    # Each hyperplane carries three masks: forbidden zero input coordinates,
    # zero state coordinates, and zero output coordinates.
    groups = {}
    families = ((tuple(tuple(int(i == j) for i in range(w)) for j in range(w)), 0),
                (injection, 1), (returned_output, 2))
    for forms, family in families:
        for index, form in enumerate(forms):
            key = canonical(form)
            masks = groups.setdefault(key, [0, 0, 0])
            masks[family] |= 1 << index
    forms, tags = tuple(groups), tuple(groups.values())
    strata = Counter()

    def add(masks, count, rank):
        if count and not masks[0]:
            strata[15 ^ masks[1], masks[2], rank] += count

    if w == 2:
        # Each distinct nonzero linear form vanishes at one projective point.
        for masks in tags:
            add(masks, 1, 1)
        add((0, 0, 0), Q+1-len(forms), 0)
    else:
        incident = {}
        for i, a in enumerate(forms):
            for j in range(i):
                b = forms[j]
                point = canonical((MUL[a[1]][b[2]] ^ MUL[a[2]][b[1]],
                                   MUL[a[2]][b[0]] ^ MUL[a[0]][b[2]],
                                   MUL[a[0]][b[1]] ^ MUL[a[1]][b[0]]))
                incident[point] = incident.get(point, 0) | (1 << i) | (1 << j)
        line_intersections = [0]*len(forms)
        for incidence in incident.values():
            masks = [0, 0, 0]
            while incidence:
                bit = incidence & -incidence
                index = bit.bit_length()-1
                line_intersections[index] += 1
                masks = [x | y for x, y in zip(masks, tags[index])]
                incidence ^= bit
            add(masks, 1, 2)
        single_count = 0
        for count, masks in zip(line_intersections, tags):
            only = Q+1-count
            assert only >= 0
            add(masks, only, 1)
            single_count += only
        add((0, 0, 0), Q*Q+Q+1-len(incident)-single_count, 0)
    assert sum(strata.values()) == (Q-1)**(w-1)
    assert all(h and zeros != 65535 for h, zeros, _ in strata)
    return dict(strata)


def projective_histogram(injection, returned_output):
    histogram = Counter()
    for (h, zeros, _), count in projective_strata(injection, returned_output).items():
        histogram[h, 16-zeros.bit_count()] += count
    return dict(histogram)


@lru_cache(maxsize=None)
def return_histogram(phase, inputs):
    return projective_histogram(*packet_forms(phase, inputs))


def joint_return_moment(phase, inputs, state_support, y):
    """Exact conditional E[y^wtY; nextstate=0] for fixed incoming support.

    Integer/Fraction y gives a Fraction.  A floating y is an exploratory
    evaluation only.  The exact histogram remains integer-valued.
    """
    if state_support not in range(1, 16) or not 0 < y <= 1:
        raise ValueError("require nonempty incoming state support and 0<y<=1")
    if isinstance(y, int):
        y = Fraction(y)
    numerator = sum(count*y**weight for (h, weight), count in return_histogram(phase, tuple(inputs)).items()
                    if h == state_support)
    return numerator/(Q-1)**(len(inputs)+state_support.bit_count()-1)


def uniform_return_bound(weight, y, *, phase=None, average_inputs=False):
    """Uniform exact-moment upper bound over incoming states and placements.

    Default maximizes over every input placement.  average_inputs=True
    first averages uniform placements for each FIXED incoming support H,
    then maximizes over H.  This requires the placement to be sampled fresh
    independently of the incoming state, conditional on its packet weight.
    No averaging over phases is performed; phase=None maximizes them.
    """
    if not 0 < y <= 1:
        raise ValueError("require 0<y<=1")
    if isinstance(y, int):
        y = Fraction(y)
    powers = [y**v for v in range(17)]
    return max(sum(count*powers[v] for v, count in profile)
               for profile in uniform_return_profiles(weight, phase, average_inputs))


@lru_cache(maxsize=None)
def uniform_return_profiles(weight, phase=None, average_inputs=False):
    """Distinct exact rational output polynomials used by the uniform bound.

    Cache the combinatorics, not a particular numeric mark.  Repeated proof
    screens therefore evaluate only the distinct short coefficient vectors.
    """
    if weight not in (1, 2, 3):
        raise ValueError("supported input weights are1,2,3")
    if phase is not None and phase not in range(4):
        raise ValueError("invalid phase")
    phases = range(4) if phase is None else (phase,)
    profiles = set()
    for current in phases:
        if average_inputs:
            totals = [Counter() for _ in range(16)]
            for inputs in combinations(range(16), weight):
                for (h, v), count in return_histogram(current, inputs).items():
                    totals[h][v] += count
            for h in range(1, 16):
                denominator = comb(16, weight)*(Q-1)**(weight+h.bit_count()-1)
                profiles.add(tuple((v, Fraction(count, denominator)) for v, count in sorted(totals[h].items())))
        else:
            for inputs in combinations(range(16), weight):
                by_state = [Counter() for _ in range(16)]
                for (h, v), count in return_histogram(current, inputs).items():
                    by_state[h][v] += count
                for h in range(1, 16):
                    denominator = (Q-1)**(weight+h.bit_count()-1)
                    profiles.add(tuple((v, Fraction(count, denominator)) for v, count in sorted(by_state[h].items())))
    return tuple(sorted(profiles))


def simple_placement_bound(weight, y):
    """Proved coarse fixed-word bound after a fresh uniform input placement."""
    if not 0 < y <= 1 or weight not in (1, 2, 3):
        raise ValueError("require weight1..3 and 0<y<=1")
    if isinstance(y, int):
        y = Fraction(y)
    if weight == 1:
        return y**16/255**4
    if weight == 2:
        return max(y**13, Fraction(17, 4)*y**16)/255**4
    return y**13/255**4


@lru_cache(maxsize=None)
def rank_return_terms(phase, inputs):
    """Integer coefficients (pivot exponent, selected-zero count, count).

    SUM over incoming witness supports H, not their maximum.  Return and
    selected output-zero equations are imposed jointly BEFORE any current
    packet diagonals are used elsewhere.  No projective-value multiplicity
    is included: the arrangement is used solely as an exact feasibility
    oracle for the selected homogeneous equation systems.
    """
    injection, output = packet_forms(phase, inputs)
    feasible = set()
    for h, zeros, _ in projective_strata(injection, output):
        subset = zeros
        while True:
            feasible.add((h, subset))
            if subset == 0:
                break
            subset = (subset-1) & zeros
    terms = Counter()
    w = len(inputs)
    for h, zeros in feasible:
        rows = [row for j, row in enumerate(injection) if not h & (1 << j)]
        rows.extend(row for j, row in enumerate(output) if zeros & (1 << j))
        rank = matrix_rank(rows, w)
        assert rank < w and zeros.bit_count() <= 3
        terms[h.bit_count()+rank, zeros.bit_count()] += 1
    return tuple((exponent, zeros, count) for (exponent, zeros), count in sorted(terms.items()))


@lru_cache(maxsize=None)
def rank_return_profiles(weight, phase=None, average_inputs=False):
    if weight not in (1, 2, 3) or (phase is not None and phase not in range(4)):
        raise ValueError("require input weight1..3 and a valid phase")
    profiles = set()
    for current in range(4) if phase is None else (phase,):
        total = Counter()
        for inputs in combinations(range(16), weight):
            terms = rank_return_terms(current, inputs)
            if average_inputs:
                for exponent, zeros, count in terms:
                    total[exponent, zeros] += count
            else:
                profiles.add(terms)
        if average_inputs:
            profiles.add(tuple((exponent, zeros, Fraction(count, comb(16, weight)))
                               for (exponent, zeros), count in sorted(total.items())))
    return tuple(sorted(profiles))


def rank_return_moment(weight, y, tau, *, phase=None, average_inputs=False):
    """Selected-zero upper action for existential alive-rank potential.

    Maximize over phase/input placement (or average fresh uniform placements)
    after SUMMING the H branches.  Selected-zero marks are
    y**(16-z)*(1-y)**z.  This is not the exact fixed-word return moment.
    """
    if not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError("require 0<y<=1 and1<=tau<=255")
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    return max(sum(count*y**(16-z)*(1-y)**z/tau**exponent
                   for exponent, z, count in profile)
               for profile in rank_return_profiles(weight, phase, average_inputs))


@lru_cache(maxsize=None)
def coarse_rank_return_profiles(weight, phase=None, average_inputs=False):
    """Return-first exact-support cover with explicit post-return rank loss.

    An entry (charge, minimum_output_weight, count) represents
    count*tau^(-charge)*y^(minimum_output_weight).  Unlike the selected-zero
    polynomial this covers exact output supports by their equation spans.
    """
    if weight not in (1, 2, 3) or (phase is not None and phase not in range(4)):
        raise ValueError("require input weight1..3 and a valid phase")
    profiles = set()
    for current in range(4) if phase is None else (phase,):
        total = Counter()
        for inputs in combinations(range(16), weight):
            injection, _ = packet_forms(current, inputs)
            hist = return_histogram(current, inputs)
            terms = Counter()
            for h in {h for h, _ in hist}:
                rows = [row for j, row in enumerate(injection) if not h & (1 << j)]
                r = matrix_rank(rows, weight)
                dimension = weight-r
                minimum = min(v for mask, v in hist if mask == h)
                for k in range(dimension):
                    terms[h.bit_count()+r+k, minimum] += comb(16, k)
            if average_inputs:
                total.update(terms)
            else:
                profiles.add(tuple((charge, minimum, count)
                                   for (charge, minimum), count in sorted(terms.items())))
        if average_inputs:
            profiles.add(tuple((charge, minimum, Fraction(count, comb(16, weight)))
                               for (charge, minimum), count in sorted(total.items())))
    return tuple(sorted(profiles))


def coarse_rank_return_moment(weight, y, tau, *, phase=None, average_inputs=False):
    """Coarse exact-support rank continuation, SUM over incoming supports H.

    Each H contributes tau^(-h-r)*y^m*sum(k<weight-r) C(16,k)*tau^(-k).
    This includes output-support branching; bare y^m alone would omit it.
    """
    if not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError("require 0<y<=1 and1<=tau<=255")
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    return max(sum(count*y**minimum/tau**charge for charge, minimum, count in profile)
               for profile in coarse_rank_return_profiles(weight, phase, average_inputs))


def summarize(weight):
    supports = Counter()
    min_exponent = 99
    minimum_weight = 16
    worst_probability = Fraction(0)
    witnesses = []
    for phase in range(4):
        for inputs in combinations(range(16), weight):
            injection, _ = packet_forms(phase, inputs)
            hist = return_histogram(phase, inputs)
            for h in sorted({h for h, _ in hist}):
                outside = [row for j, row in enumerate(injection) if not h & (1 << j)]
                exponent = h.bit_count()+matrix_rank(outside, weight)
                minimum_weight = min(minimum_weight, min(v for (mask, v) in hist if mask == h))
                min_exponent = min(min_exponent, exponent)
                supports[phase, h.bit_count(), exponent] += 1
                probability = joint_return_moment(phase, inputs, h, Fraction(1))
                assert probability <= Fraction(1, 255**exponent)
                if probability > worst_probability:
                    worst_probability = probability
                    witnesses = [(phase, inputs, h)]
                elif probability == worst_probability:
                    witnesses.append((phase, inputs, h))
    return dict(weight=weight, minimum_pivot_exponent=min_exponent,
                minimum_return_output_weight=minimum_weight,
                worst_probability=str(worst_probability), witnesses= witnesses[:8],
                number_of_worst_witnesses=len(witnesses),
                feasible_support_counts=[dict(phase=p, state_weight=h, exponent=e, count=n)
                                         for (p, h, e), n in sorted(supports.items())])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weights", type=int, nargs="+", default=(2, 3), choices=(1, 2, 3))
    args = parser.parse_args()
    for weight in args.weights:
        print(json.dumps(summarize(weight)), flush=True)
