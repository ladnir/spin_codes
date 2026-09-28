"""Outward shape-resolved eleven-coordinate IMT bounds for independent rows.

The state coordinates retain the existing mature-tail semantics. Shapes
are averaged only after all their coupled coefficients are constructed.
Above the configured cutoff, a valid coarse universal envelope is used.
This module supplies local operators, not a complete distance certificate.
"""
import argparse
from collections import Counter
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial, prod
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
from flint import arb, arb_mat, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from group_rank_one_verify import up
from occupancy_model import local_data
from occupancy_memory import (prepare as prepare_memory, epoch_operators, coarse_epoch,
                              rounded, Z, F, M, C, U)
from occupancy_fresh_moment import fresh_census
from occupancy_window_average import prepare_inputs, averages
from occupancy_multi_average import multiplicities
from fresh_collision import probability as fresh_collision_probability
from zero_moment import census as zero_census
from mature_tail import census as tail_census, L48, L56, TAIL_TERMINAL
from pair_tail import census as pair_census
from mixing_rounds import transform
import full_feedback_refinement
import window_histogram
import cancellation_joint


def rational(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def weights_of(counts):
    return tuple(b for b, n in enumerate(counts, 1) for _ in range(n))


def prepare_shared(args):
    """Exact censuses shared across tilts; no binary64 cache is consumed."""
    cut = getattr(args, 'cut', 8)
    full = min(cut, getattr(args, 'full_feedback', 6))
    window = min(cut, getattr(args, 'window_histogram', 8))
    joint = getattr(args, 'joint_cancellation', True)
    prepared = prepare_memory(local_data(4))
    inputs = prepare_inputs()
    data = SimpleNamespace(prepared=prepared, inputs=inputs,
        fresh=fresh_census(prepared), tails=tail_census(inputs, prepared),
        pairs=pair_census(inputs, prepared), zeros=zero_census(),
        full_max=full, window_max=window)
    data.feedback = full_feedback_refinement.census(full) if full >= 2 else {}
    if data.feedback:
        full_feedback_refinement.check_pairs(data.feedback)
    data.histograms = window_histogram.census(prepared) if window else None
    data.cancellations = cancellation_joint.census(check_old=False) if joint else None
    if joint:
        cancellation_joint.check_fresh(data.cancellations, prepared)
        if full >= 3:
            cancellation_joint.check_feedback(data.cancellations, data.feedback)
        for shape, (rows, _, denominator) in data.cancellations[0].items():
            if len(shape) >= 2:
                old, old_den = data.zeros[len(shape)][shape]
                assert all(sum(row[w] for row in rows.values())*old_den == old[w]*denominator
                           for w in range(129))
        print('Shape-inner joint marginals match independently enumerated output counts', flush=True)
    density_maximum = getattr(args, 'feedback_density', 0)
    data.feedback_density = None
    if density_maximum:
        if not 1 <= density_maximum <= min(6,cut):
            raise ValueError('feedback-density census requires maximum 1..min(6,cut)')
        from feedback_density import build as build_density
        tilts = getattr(args, 'tilts', None) or [args.tilt]
        data.feedback_density = build_density(density_maximum,tilts,getattr(args,'precision',192),2)
        print('Feedback-density provenance:',data.feedback_density['provenance'],flush=True)
    return data


def coarse_shape(spectrum, weights, powers, pair_bound):
    """The existing coarse nine-coordinate proof for one fixed shape."""
    j, weight = len(weights), sum(weights)
    m = (1 << 19)-1
    atom = min(arb(1)/((33-j)*max(comb(4, b) for b in weights)),
               rational(pair_bound)*32*31/((34-j)*(33-j)))
    f = powers[max(0, 48-weight)]
    t = arb_mat(9, 9)
    t[Z,Z], t[Z,M], t[Z,C] = powers[weight]*atom, powers[weight], powers[weight]*atom
    fresh = f*min(atom, arb(1)/32)
    t[F,Z], t[F,M], t[F,C] = fresh/2+f/(2*m), f/2, fresh/2
    t[M,Z], t[M,M] = f/(2*m), f/2
    t[C,Z], t[C,C] = f/2, f/2
    for source in (F, M):
        for k, v in enumerate(sorted(spectrum)):
            t[source,U+k] = f*spectrum[v]/(2*m)
    for i, v in enumerate(sorted(spectrum)):
        peak = powers[max(0, v-weight)]
        t[U+i,Z] = peak/(2*spectrum[v])+peak/(2*m)
        t[U+i,M], t[U+i,C] = peak/2, peak/(2*spectrum[v])
        for k, w in enumerate(sorted(spectrum)):
            t[U+i,U+k] = peak*spectrum[w]/(2*m)
    return rounded(t)


def refine_nine(t, weights, data, window, powers, factors):
    """Per-shape counterparts of the existing nine-coordinate refinements."""
    t = t*1
    j, weight = len(weights), sum(weights)
    spectrum, distributions, fresh_moments = data.fresh
    levels, m = sorted(spectrum), (1 << 19)-1
    def put(i, k, value):
        t[i,k] = min(t[i,k], up(value))
    def average(hist, shift=0):
        return up(sum((count*powers[max(0, v-shift)] for v, count in hist.items()), arb(0))/sum(hist.values()))
    if j == 1:
        moment = max(average(hist) for (a,b), hist in fresh_moments.items() if b == weights[0])
    else:
        moment = max(average(hist, weight) for hist in distributions.values())
    if j:
        put(F, M, moment/2)
    for k, v in enumerate(levels):
        put(F, U+k, moment*spectrum[v]/(2*m))
    if j == 1:
        mass, density, zero, denominator = window[weights[0]]
        mass, density, zero = (arb(x)/denominator for x in (mass, density, zero))
        put(M, M, mass/2); put(M, Z, mass/(2*m))
        put(C, Z, zero/2); put(C, C, density/2)
        for k, v in enumerate(levels):
            put(M, U+k, mass*spectrum[v]/(2*m))
            put(U+k, C, density/(2*spectrum[v]))
    if j >= 2:
        moments = {v:min(arb(1), up(powers[v]*prod(factors[v][b-1] for b in weights))) for v in levels}
        arbitrary = max(moments.values())
        fresh = max(up(sum((n*moments[v] for v,n in hist.items()), arb(0))/sum(hist.values()))
                    for hist in distributions.values())
        put(F, M, fresh/2); put(M, M, arbitrary/2); put(M, Z, arbitrary/(2*m))
        for k, v in enumerate(levels):
            put(F, U+k, fresh*spectrum[v]/(2*m))
            put(M, U+k, arbitrary*spectrum[v]/(2*m))
            put(U+k, M, moments[v]/2)
            put(U+k, Z, powers[abs(v-weight)]/(2*spectrum[v])+moments[v]/(2*m))
            for ell, w in enumerate(levels):
                put(U+k, U+ell, moments[v]*spectrum[w]/(2*m))
    if j in (2,3):
        tables = {1:{}, 2:{}, **data.prepared[0][2]}
        cancellation = max(fresh_collision_probability(tables, a, weights) for a in range(1,5))
        put(F, Z, powers[max(0,48-weight)]*(rational(cancellation)+arb(1)/m)/2)
        hist, denominator = data.zeros[j][weights]
        put(C, Z, sum((n*p for n,p in zip(hist,powers)),arb(0))/(2*denominator))
    return rounded(t)


def lift_shape(old, weights, data, powers, factors):
    """Mature-tail lift with the (a,b,c) split computed for this shape.

    L48 and L56 remain upper masses of subsets of M, not extra mass.
    All coefficients initially describe one transvection; transform()
    changes them to two only after this lift.
    """
    j, weight = len(weights), sum(weights)
    spectrum = data.inputs[3]
    levels, m = sorted(spectrum), (1 << 19)-1
    t = arb_mat([[old[i,k] if i < 9 and k < 9 else arb(0) for k in range(11)] for i in range(11)])
    for source in [Z,F]+list(range(U,U+5)):
        t[source,L48] = old[source,M]
        t[source,L56] = old[source,M]
    moments = {v:min(arb(1), up(powers[v]*prod(factors[v][b-1] for b in weights))) for v in levels}
    a = max(moments[v] for v in levels if v >= 64)
    b = max(arb(0), up(moments[56]-a))
    c = max(arb(0), up(moments[48]-a-b))
    for source, factor in ((M,a),(L56,b),(L48,c)):
        t[source,M] = up(factor/2)
        if j:
            t[source,Z] = up(factor/(2*m))
        for k, v in enumerate(levels):
            t[source,U+k] = up(factor*spectrum[v]/(2*m))
    if not j:
        t[L48,L48] = powers[48]/2
        t[L56,L56] = powers[56]/2
        # This difference need not upper-bound the difference of the exact
        # exponentials separately. The paired coefficients sum to at least
        # powers[48] on L48, while the L56-only coefficient bounds weight56.
        # Do not tighten those two coefficients independently afterward.
        t[L48,L56] = up((powers[48]-powers[56])/2)
    elif j == 1:
        column = weights[0]
        for target, cutoff in ((L48,48),(L56,56)):
            arbitrary, fresh, uniform = data.tails[column][cutoff][:3]
            for source in [F,M]+list(range(U,U+5)):
                p = arbitrary if source == M else fresh if source == F else uniform[levels[source-U]]
                v = 48 if source in (F,M) else levels[source-U]
                t[source,target] = min(old[source,M], up(powers[max(0,v-column)]*rational(p)/2))
    else:
        # This is the pair-tail refinement in the original lift, now with
        # one fixed shape rather than a maximum over incompatible shapes.
        choices = set(combinations(weights,2))
        for target, cutoff in ((L48,48),(L56,56)):
            probability = min(Q(1), min(data.pairs[p][cutoff][0] for p in choices)
                              *Q(32*31,(34-j)*(33-j)))
            zero_probability = data.pairs[weights][cutoff][1] if j == 2 else probability
            for source in [Z,F,M]+list(range(U,U+5)):
                if source == Z:
                    value = powers[weight]*rational(zero_probability)
                else:
                    v = 48 if source in (F,M) else levels[source-U]
                    value = powers[max(0,v-weight)]*rational(probability)/2
                t[source,target] = min(old[source,M], up(value))
            t[L48,target] = t[L56,target] = arb(0)
    return rounded(t)


def coarse_lift(old):
    """Embed a valid mature-mass bound by giving each tail its whole mass."""
    t = arb_mat([[old[i,k] if i < 9 and k < 9 else arb(0) for k in range(11)] for i in range(11)])
    for source in range(9):
        t[source,L48] = old[source,M]
        t[source,L56] = old[source,M]
    return rounded(t)


def transform_one(t, spectrum, tilt, occupied):
    # transform distinguishes empty from nonempty occupancy, not its value.
    return transform([arb_mat(11,11), t] if occupied else [t], spectrum, tilt, 2)[-1]


def refine_one(function, t, j, *args):
    """Use an existing refinement with a singleton shape census.

    Only slot j is visited by the supplied one-shape data. The unused zero
    matrices are indexing placeholders, never transfer operators.
    """
    base = [arb_mat(11,11) for _ in range(j)]+[t]
    # The legacy functions announce each refinement. Here they run once
    # per shape; retain one occupancy summary instead of hundreds of lines.
    with redirect_stdout(StringIO()):
        result = function(base, *args)[j]
    return result


class ShapeInner:
    def __init__(self, tilt, precision, cut, shapes, fallback, shared):
        self.tilt, self.precision, self.cut = tilt, precision, cut
        self.shapes, self.fallback, self.shared = shapes, fallback, shared
        self.terminal = TAIL_TERMINAL
        self.float_shapes = {}
        for j, values in shapes.items():
            counts = list(values)
            matrices = np.array([[[float(t[i,k]) for k in range(11)] for i in range(11)] for t in values.values()])
            coefficients = np.array([factorial(j)//prod(factorial(n) for n in c) for c in counts], dtype=float)
            self.float_shapes[j] = np.array(counts, dtype=int), matrices, coefficients
        self.float_fallback = [np.array([[float(t[i,k]) for k in range(11)] for i in range(11)]) for t in fallback]

    def mix(self, theta):
        """Outward local operators for exact rational nonzero weight law."""
        theta = tuple(Q(x) for x in theta)
        if len(theta) != 4 or min(theta) < 0 or sum(theta) != 1:
            raise ValueError('four nonnegative rational probabilities summing to one required')
        ctx.prec = self.precision
        result = []
        for j in range(33):
            if j > self.cut:
                result.append(self.fallback[j]*1)
                continue
            matrix = arb_mat(11,11)
            total = Q(0)
            for counts, t in self.shapes[j].items():
                mass = Q(factorial(j),prod(factorial(n) for n in counts))*prod(p**n for p,n in zip(theta,counts))
                total += mass
                if mass:
                    matrix += rational(mass)*t
            assert total == 1
            result.append(rounded(matrix))
        return result

    def float_mix(self, theta):
        """Binary64 proposals only; outward callers must use mix()."""
        theta = np.asarray(theta, dtype=float)
        if theta.shape != (4,) or np.min(theta) < 0 or abs(theta.sum()-1) > 1e-12:
            raise ValueError('invalid proposal probabilities')
        result = []
        for j in range(33):
            if j > self.cut:
                result.append(self.float_fallback[j])
            else:
                counts, matrices, coefficient = self.float_shapes[j]
                mass = coefficient*np.prod(theta[None,:]**counts,axis=1)
                assert abs(mass.sum()-1) < 1e-12
                result.append(np.einsum('s,sik->ik',mass,matrices,optimize=False))
        return np.array(result)


def build(args, shared=None):
    """Build for one tilt; reuse prepare_shared(args) across different tilts."""
    tilt, precision = str(args.tilt), getattr(args, 'precision', 192)
    cut, penalty = getattr(args, 'cut', 8), Q(getattr(args, 'penalty', 1))
    if not 4 <= cut <= 10 or penalty != 1 or precision < 128 or Q(tilt) <= 0:
        raise ValueError('requires cut4..10, penalty1, precision>=128, positive tilt')
    shared = prepare_shared(args) if shared is None else shared
    ctx.prec = precision
    spectrum = shared.inputs[3]
    powers = [up((-arb(tilt)*w).exp()) for w in range(145)]
    factors = {v:[(arb(128-v)*(-arb(tilt)*b).exp()+arb(v)*(arb(tilt)*b).exp())/128
                  for b in range(1,5)] for v in spectrum}
    raw = epoch_operators(shared.prepared, tilt, detailed=True)
    window = averages(shared.inputs, tilt)
    histogram_moments = (window_histogram.moments(shared.histograms, tilt, min(cut,shared.window_max))
                         if shared.histograms is not None else {})
    pair_bound = max(Q(row[2],row[0]) for row in shared.prepared[0][1][1].values())
    shapes = {}
    for j in range(cut+1):
        shapes[j] = {}
        for counts in multiplicities(j):
            weights = weights_of(counts)
            old = raw[0] if j == 0 else raw[j][weights] if j <= 4 else coarse_shape(spectrum, weights, powers, pair_bound)
            old = refine_nine(old, weights, shared, window, powers, factors)
            t = transform_one(lift_shape(old, weights, shared, powers, factors), spectrum, tilt, bool(j))
            if j and counts in histogram_moments:
                t = refine_one(window_histogram.refine, t, j, shared.histograms,
                               {counts:histogram_moments[counts]}, 1, 2)
            if j >= 2 and weights in shared.feedback:
                t = refine_one(full_feedback_refinement.refine, t, j, {weights:shared.feedback[weights]}, spectrum, tilt, 1, 2)
            if j and shared.cancellations is not None and weights in shared.cancellations[0]:
                subset = ({weights:shared.cancellations[0][weights]}, spectrum)
                t = refine_one(cancellation_joint.refine, t, j, subset, tilt, 1, 2)
            t = rounded(t)
            if getattr(args, 'class_tail', False):
                from class_tail import refine as refine_class_tail
                t = refine_class_tail(t, weights, shared, powers, 2)
            if getattr(args, 'column_density', False):
                from column_density import refine as refine_column_density
                t = refine_column_density(t, weights, window, tilt, spectrum=spectrum, rounds=2)
            if getattr(shared, 'feedback_density', None) is not None:
                density = shared.feedback_density['bounds'][tilt].get(weights)
                if density is not None:
                    assert all(t[source,C] == 0 for source in (M,L48,L56))
                    t[C,C] = min(t[C,C],density['density'])
                    for k,v in enumerate(sorted(spectrum)):
                        t[U+k,C] = min(t[U+k,C],density['uniform'][v])
            assert all(t[i,k] >= 0 for i in range(11) for k in range(11))
            shapes[j][counts] = t
        print('Outward shape-inner tilt/local occupancy/shapes:', tilt, j, len(shapes[j]), flush=True)
    fallback = [arb_mat(11,11) for _ in range(cut+1)]
    for j in range(cut+1,33):
        old = coarse_epoch(spectrum, j, tilt, pair_bound=pair_bound)
        fallback.append(transform_one(coarse_lift(old), spectrum, tilt, True))
    result = ShapeInner(tilt, precision, cut, shapes, fallback, shared)
    self_test(result)
    return result


def self_test(model):
    """Exact mixture identities and direct state checks through three windows."""
    checks = 0
    for b in range(4):
        theta = tuple(Q(int(i == b)) for i in range(4))
        exact, proposal = model.mix(theta), model.float_mix(theta)
        for j in range(model.cut+1):
            counts = tuple(j*int(i == b) for i in range(4))
            target = model.shapes[j][counts]
            for i in range(11):
                for k in range(11):
                    assert exact[j][i,k] >= target[i,k]
                    assert abs(float(exact[j][i,k])-proposal[j,i,k]) < 1e-12
                    checks += 1
    assert all(model.terminal[i] == 0 for i in (C,L48,L56))
    print('Shape-inner atomic mixture/coordinate checks:', checks, 'passed', flush=True)
    direct_state_test(model)


def direct_state_test(model):
    """Independent exact integer input census and higher-precision moments."""
    low, high, words, spectrum = model.shared.inputs
    expansion = np.bitwise_count(low)+np.bitwise_count(high)
    states = sorted(set([1,17]+[int(np.flatnonzero(expansion == v)[0]) for v in sorted(spectrum)]))
    cases = {():[(0,0,0)]}
    cases.update({(b,):entries for b,entries in words.items()})
    for shape in ((1,1),(1,4),(2,2),(4,4)):
        a,b = shape
        values = []
        for sa,la,ha in words[a]:
            wa = ((la | (ha << 64)).bit_length()-1)//4
            for sb,lb,hb in words[b]:
                wb = ((lb | (hb << 64)).bit_length()-1)//4
                if wa != wb:
                    values.append((sa ^ sb, la ^ lb, ha ^ hb))
        assert len(values) == 32*31*comb(4,a)*comb(4,b)
        cases[shape] = values
    by_window = {b:[[] for _ in range(32)] for b in range(1,5)}
    for b, entries in words.items():
        for syndrome, lo, hi in entries:
            window = ((lo | (hi << 64)).bit_length()-1)//4
            by_window[b][window].append((syndrome,lo,hi))
    for shape in ((4,4,4),(1,4,4)):
        values = []
        for windows in combinations(range(32),3):
            for assignment in set(permutations(shape)):
                for chosen in product(*(by_window[b][w] for b,w in zip(assignment,windows))):
                    values.append(tuple(chosen[0][i] ^ chosen[1][i] ^ chosen[2][i] for i in range(3)))
        multiplicity = factorial(3)//prod(factorial(n) for n in Counter(shape).values())
        assert len(values) == comb(32,3)*multiplicity*prod(comb(4,b) for b in shape)
        cases[shape] = values
    old_precision = ctx.prec
    ctx.prec = model.precision+128
    powers = [(-arb(model.tilt)*w).exp() for w in range(129)]
    m, checks = (1 << 19)-1, 0
    try:
        for shape, entries in cases.items():
            counts = tuple(shape.count(b) for b in range(1,5))
            t = model.shapes[len(shape)][counts]
            syndromes = np.array([s for s,_,_ in entries],dtype=np.uint64)
            lows = np.array([lo for _,lo,_ in entries],dtype=np.uint64)
            highs = np.array([hi for _,_,hi in entries],dtype=np.uint64)
            for state in states:
                weights = np.bitwise_count(np.uint64(low[state]) ^ lows)+np.bitwise_count(np.uint64(high[state]) ^ highs)
                targets = syndromes ^ np.uint64(state)
                def moment(selected):
                    hist = np.bincount(weights[selected],minlength=129)
                    return sum((int(n)*p for n,p in zip(hist,powers)),arb(0))/len(entries)
                mass = moment(np.ones(len(entries),dtype=bool))
                actual = [arb(0) for _ in range(11)]
                actual[M] = moment(targets != 0)/4
                actual[Z] = moment(targets == 0)/4+3*moment(syndromes != 0)/(4*m)
                for target,cutoff in ((L48,48),(L56,56)):
                    actual[target] = moment((targets != 0)&(expansion[targets] <= cutoff))/4
                for k,v in enumerate(sorted(spectrum)):
                    actual[U+k] = 3*mass*spectrum[v]/(4*m)
                histograms = {}
                for target,w in zip(targets,weights):
                    if target:
                        histograms.setdefault(int(target),Counter())[int(w)] += 1
                actual[C] = max((sum((n*powers[w] for w,n in hist.items()),arb(0))/(4*len(entries))
                                 for hist in histograms.values()),default=arb(0))
                incoming = [0]*11
                incoming[M] = incoming[C] = 1
                incoming[L48] = int(expansion[state] <= 48)
                incoming[L56] = int(expansion[state] <= 56)
                for k, value in enumerate(actual):
                    upper = sum((t[i,k]*incoming[i] for i in range(11)),arb(0))
                    assert value <= upper or (value-upper).contains(0), (shape,state,k,value,upper)
                    checks += 1
    finally:
        ctx.prec = old_precision
    print('Shape-inner direct empty/single/two/three-window state inequalities:', checks, 'passed', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt', default='.032')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--cut', type=int, default=8)
    parser.add_argument('--full-feedback', type=int, default=6)
    parser.add_argument('--window-histogram', type=int, default=8)
    parser.add_argument('--no-joint-cancellation', action='store_false', dest='joint_cancellation')
    build(parser.parse_args())
