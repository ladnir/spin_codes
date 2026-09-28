"""Weighted fresh-state envelope: diagnostic only, not a full certificate.

For fresh feedback distributions p_a, use generators p_a / kappa_a,
kappa_a >= 1. Existing fresh-source coefficients remain upper bounds.
Activation into this cone costs kappa_b; selected outgoing coefficients
can instead divide each incoming-shape bound by its own kappa_a.
"""
import argparse
from fractions import Fraction
from math import comb, log, prod

import numpy as np
from flint import arb, ctx
from scipy.optimize import minimize_scalar

from group_rank_one_verify import up
from occupancy_memory import Z, F, M, C, U, prepare
from occupancy_model import local_data
from occupancy_fresh_moment import fresh_census
from fresh_collision import probability
from occupancy_multi_average import multiplicities
from occupancy_sensitivity import float_placement
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment, log_binomial_mass
from mature_tail import TAIL_TERMINAL


def coefficients(prepared, census, tilt, penalty, gauge, rounds=2, degree=4):
    """Arb upper coefficients, with the old cone recovered at kappa=1."""
    assert len(gauge) == 4 and min(Fraction(x) for x in gauge) >= 1
    kappa = {a: arb(str(x)) for a, x in enumerate(gauge, 1)}
    spectrum, levels, moments = census
    single, pair, zeros = prepared[0]
    assert all(row[1] == 0 for row in pair[1].values())
    cancel = single[4]
    m = (1 << 19) - 1
    lazy = arb(2) ** (-rounds)
    refresh = 1 - lazy
    powers = [(-arb(tilt) * w).exp() for w in range(129)]
    rho = arb(penalty)

    def mean(hist):
        return sum((c * powers[w] for w, c in hist.items()), arb(0)) / sum(hist.values())

    result = {}
    # Output weighting preserves a pointwise cap by the same fresh generator
    # during an empty lazy step. Do not replace that self-loop by an average.
    empty = max(up(mean(hist) / kappa[a]) for a, hist in levels.items())
    for k, w in enumerate(sorted(spectrum)):
        result[0, F, U + k] = up(refresh * empty * spectrum[w] / m)
    result[1, Z, F] = max(up(powers[b] * rho ** int(b == 4) * kappa[b]) for b in range(1, 5))

    active_moment = max(up(mean(hist) * rho ** int(b == 4) / kappa[a])
                        for (a, b), hist in moments.items())
    result[1, F, M] = up(lazy * active_moment)
    for k, w in enumerate(sorted(spectrum)):
        result[1, F, U + k] = up(refresh * active_moment * spectrum[w] / m)

    zero_terms = []
    density_terms = []
    for a in range(1, 5):
        for b in range(1, 5):
            choices = 32 * comb(4, b)
            # Disjoint fresh supports: p_a can cancel this input only if a=b.
            cancellation = (sum((count * powers[w] for row in cancel[(b,)].values()
                                 for w, count in row.items()), arb(0)) / choices**2
                            if a == b else arb(0))
            scale = rho ** int(b == 4) / kappa[a]
            zero_terms.append(up(scale * (lazy * cancellation + refresh * mean(moments[a, b]) / m)))
            peak = prepared[1][tuple(sorted(((a,), (b,))))][1]
            density_terms.append(up(scale * lazy * powers[48-b] * arb(peak.numerator) / peak.denominator))
    result[1, F, Z] = max(zero_terms)
    result[1, F, C] = max(density_terms)

    # For more inputs retain either the shift bound or the proven
    # without-replacement chord bound, before maximizing over fresh types.
    tables = {1: {}, 2: {}, **zeros}
    for j in range(2, degree + 1):
        masses = []
        cancellations = []
        for shape in multiplicities(j):
            weights = tuple(b for b, count in enumerate(shape, 1) for _ in range(count))
            total = sum(weights)
            scale = rho ** shape[3]
            by_weight = {v: min(arb(1), up(powers[v] * prod(
                ((arb(128-v) * powers[b] + arb(v) / powers[b]) / 128) ** count
                for b, count in enumerate(shape, 1)))) for v in spectrum}
            for a, hist in levels.items():
                shift = sum((c * powers[max(0, v-total)] for v, c in hist.items()), arb(0)) / sum(hist.values())
                chord = sum((c * by_weight[v] for v, c in hist.items()), arb(0)) / sum(hist.values())
                moment = min(up(shift), up(chord))
                masses.append(up(scale * moment / kappa[a]))
                if j in (2, 3):
                    prob = probability(tables, a, weights)
                    cancellation = powers[max(0, 48-total)] * arb(prob.numerator) / prob.denominator
                    cancellations.append(up(scale / kappa[a] * (lazy*cancellation + refresh*moment/m)))
        mass = max(masses)
        result[j, F, M] = up(lazy * mass)
        for k, w in enumerate(sorted(spectrum)):
            result[j, F, U+k] = up(refresh * mass * spectrum[w] / m)
        if cancellations:
            result[j, F, Z] = max(cancellations)
    return result


def apply(base, entries):
    result = base.copy()
    for (j, source, target), value in entries.items():
        value = float(value)
        if (j, source, target) == (1, Z, F):
            result[j, source, target] = value
        else:
            result[j, source, target] = min(result[j, source, target], value)
    return result


def self_test(prepared, census):
    # Independent direct summation over each fresh atom and each single input.
    # Test the lazy zero contribution and total lazy mass, plus the refresh.
    from group_moment import maps
    images, columns, spectrum = maps()
    atoms = prepared[0][0][3]
    cases = 0
    for gauge in (('1','1','1','1'), ('1','1.2','1.5','2')):
        tilt = '.032'; penalty = '.75'
        entries = coefficients(prepared, census, tilt, penalty, gauge, degree=3)
        for a in range(1, 5):
            for b in range(1, 5):
                mass = arb(0); zero = arb(0)
                for state, count in atoms[(a,)].items():
                    for window in range(32):
                        for bits in range(1, 16):
                            if bits.bit_count() != b:
                                continue
                            feedback = 0
                            for bit in range(4):
                                if bits >> bit & 1:
                                    feedback ^= columns[4*window+bit]
                            weight = (images[state] ^ (bits << (4*window))).bit_count()
                            value = count * (-arb(tilt)*weight).exp()
                            mass += value
                            if state == feedback:
                                zero += value
                denominator = sum(atoms[(a,)].values()) * 32 * comb(4,b)
                scale = arb(penalty)**int(b==4) / arb(gauge[a-1]) / denominator
                for actual, bound in ((scale*(zero/4 + 3*mass/(4*((1<<19)-1))), entries[1,F,Z]),
                                      (scale*mass/4, entries[1,F,M])):
                    assert actual <= bound or (actual-bound).contains(0)
                assert entries[1,Z,F] >= (-arb(tilt)*b).exp()*arb(penalty)**int(b==4)*arb(gauge[b-1])
                cases += 1
    print('Fresh gauge: independent single-input checks:', cases, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', default='tmp/r2-q64-sensitivity.npz')
    parser.add_argument('--gauge', nargs=4, action='append')
    args = parser.parse_args()
    ctx.prec = 192
    with np.load(args.snapshot, allow_pickle=False) as saved:
        base = saved['base']; outer = float(saved['outer'])
        groups, support, tilt, penalty = map(str, saved['parameters'])
    groups = int(groups); support = int(support)
    assert base.shape == (min(groups,32)+1,11,11) and groups >= 4
    prepared = prepare(local_data(4)); census = fresh_census(prepared)
    self_test(prepared, census)
    gauges = args.gauge or [('1','1','1','1'), ('1','1','1','1.25'), ('1','1','1','1.5'),
                            ('1','1','1','2'), ('1','1.2','1.2','1.5'), ('1','1.5','1.5','2')]
    for gauge in [None] + gauges:
        operators = base if gauge is None else apply(base, coefficients(prepared, census, tilt, penalty, gauge))
        region = float_placement(operators, groups)
        def objective(z):
            p = 1/(1+np.exp(-z))
            return (log_power_moment(matrix_for_probabilities(region,[p]*groups),256,TAIL_TERMINAL)
                    +float(tilt)*209715-groups*log_binomial_mass(256,support,p))
        fit = minimize_scalar(objective, bounds=(-8.,16.), method='bounded')
        print('DIAGNOSTIC NOT CERTIFICATE: gauge', gauge, 'log2 score', (fit.fun+outer)/log(2),
              'p', 1/(1+np.exp(-fit.x)), flush=True)


if __name__ == '__main__':
    main()
