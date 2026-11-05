"""Positive FLOAT regional, all-profile, box, and terminal-tail screens.

The construction is the literal lift defined by hd_field128_lift_screen.
Every marker evaluation bounds every coefficient. A bank takes entrywise
minima only after normalizing the requested histogram. One state potential
then covers every alternating matching profile at a selected group weight.
Band evaluation visits every integer weight; boxes use a convex envelope
over their entire clipped histogram domain. None is a directed certificate.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from flint import arb, ctx
from scipy.optimize import minimize_scalar
from scipy.special import gammaln

from hd_field128_lift_screen import (Budget, FIELD_ORDER, GROUPS, PAIRS, REGIONS,
    OUTPUT_LENGTH, LN2, NEG, checked_q, local_logs, log1mexp2, outer_log_density,
    source_pins as sparse_source_pins)
from hd_grouped_packet_screen import OCC
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_cover_backend import histogram_vertices, profile_count, validate_box
from hd_shifted_good_pair_screen import good_geometry
from hd_shifted_pairs_screen import half_geometry


LOG_FACTORIAL = gammaln(np.arange(PAIRS+1)+1)/LN2


def full_counts(pairs, weight):
    if type(pairs) is not int or type(weight) is not int or not 1 <= pairs <= PAIRS or not 0 <= weight <= 2*pairs:
        raise ValueError('Expected integer pair count and admissible weight')
    return np.arange(max(0, weight-pairs), weight//2+1, dtype=int)


def _valid_markers(rows):
    rows = tuple(dict.fromkeys(tuple(map(float, row)) for row in rows))
    if (not rows or len(rows) > 64 or any(len(row) != 2 or not all(map(isfinite, row))
                                        or any(abs(x) > 4096 for x in row) for row in rows)):
        raise ValueError('Expected 1..64 finite marker pairs with log magnitudes <=4096')
    return rows


class PositivePacket:
    """Prepare q-dependent packet logs once for repeated marker evaluations."""

    def __init__(self, z, alpha, *, q=FIELD_ORDER, max_seconds=60., _budget=None):
        checked_q(q, construction=True)
        if not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Expected 0<z<1 and 0<alpha<=1')
        budget = _budget or Budget(max_seconds)
        budget.check()
        self.z, self.alpha, self.q = float(z), float(alpha), q
        logs = local_logs(self.z, q).reshape((-1, 4))
        halves, terms, indices = half_geometry()
        rows = [(3*i+min(2, s), s, f, n-s+(a != b))
                for i, ((a, b), row) in enumerate(zip(halves, terms)) for s, f, n in row]
        self.half_bins = np.array([row[0] for row in rows], dtype=int)
        self.half_single = np.array([row[1] for row in rows], dtype=int)
        self.half_full = np.array([row[2] for row in rows], dtype=int)
        self.half_weights = np.array([row[3] for row in rows])
        self.half_logs = self.alpha*logs[indices]
        s, f, counts, indices = good_geometry()
        self.good_s, self.good_f = s, f
        self.good_weights = np.log2(counts)
        good = logs[indices].reshape((-1, 2, 2))
        self.good_logs = self.alpha*good.reshape((-1, 4))
        composed = np.logaddexp2(good[:, None, :, 0, None]+good[None, :, None, 0, :],
                                 good[:, None, :, 1, None]+good[None, :, None, 1, :])
        self.composed_logs = self.alpha*composed.reshape((len(s), len(s), 4))
        budget.check()

    def pair_logs(self, odds, *, max_seconds=60., _budget=None):
        """GG before alpha; GD, DG, DD after their one-step alpha powers."""
        odds = _valid_markers((odds,))[0]
        budget = _budget or Budget(max_seconds)
        budget.check()
        values = self.half_weights+self.half_single*odds[0]+self.half_full*odds[1]
        halves = np.full(153*3, NEG)
        np.logaddexp2.at(halves, self.half_bins, values)
        halves = halves.reshape((153, 3))
        clean, one, dirty = halves.T
        combined = np.logaddexp2(
            dirty[:, None]+np.logaddexp2.reduce(halves, axis=1)[None, :],
            np.logaddexp2(np.logaddexp2(clean, one)[:, None]+dirty[None, :],
                         one[:, None]+one[None, :]))
        dirty_logs = np.logaddexp2.reduce(combined.reshape((-1, 1))+self.half_logs, axis=0)
        weights = self.good_weights+self.good_s*odds[0]+self.good_f*odds[1]
        good_logs = np.logaddexp2.reduce(weights[:, None]+self.good_logs, axis=0)
        composed = np.logaddexp2.reduce((weights[:, None, None]+weights[None, :, None]
                                        +self.composed_logs).reshape((-1, 4)), axis=0)
        good_matrix, dirty_matrix = good_logs.reshape((2, 2)), dirty_logs.reshape((2, 2))
        result = composed.reshape((2, 2))
        for left, right in ((good_matrix, dirty_matrix), (dirty_matrix, good_matrix),
                            (dirty_matrix, dirty_matrix)):
            result = np.logaddexp2(result, np.asarray(multiply(left, right)))
        if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
            raise ArithmeticError('Invalid positive paired matrix logarithm')
        budget.check()
        return result


def contraction_logs(regions, log_ratio):
    if not isfinite(log_ratio) or abs(log_ratio) > 4096:
        raise ValueError('Finite common log potential ratio in [-4096,4096] required')
    a, b, c, d = np.asarray(regions).reshape((-1, 4)).T
    return np.maximum(np.logaddexp2(a, b+log_ratio), np.logaddexp2(c-log_ratio, d))


def support_sum_logs(pairs, weight, regions, log_ratio, *, method='marginal', repetitions=256,
                     max_seconds=60., _budget=None):
    """All indexed support sets, with no independent-matchings assumption."""
    full = full_counts(pairs, weight)
    if np.asarray(regions).shape != (len(full), 2, 2) or repetitions < 2 or repetitions % 2:
        raise ValueError('One matrix per feasible full count and positive even repetitions required')
    budget = _budget or Budget(max_seconds)
    budget.check()
    contraction = contraction_logs(regions, log_ratio)
    if method == 'marginal':
        singles, empty = weight-2*full, pairs-weight+full
        counts = ((gammaln(pairs+1)-gammaln(empty+1)-gammaln(singles+1)
                   -gammaln(full+1))/LN2+singles)
        value = np.logaddexp2.reduce(counts+repetitions*contraction)
    elif method == 'runs':
        # Sum log terms in row chunks. Every (a,b) is included, but no
        # quadratic matrix or FFT cancellation floor is needed.
        length = 2*pairs
        runs = np.arange(weight+1)
        counts = np.full(weight+1, NEG)
        if weight in (0, length):
            counts[0] = 0.
        else:
            valid = (runs >= 1) & (runs <= min(weight, length-weight))
            c = runs[valid]
            counts[valid] = (np.log2(length/c)
                +(gammaln(weight)-gammaln(c)-gammaln(weight-c+1)
                  +gammaln(length-weight)-gammaln(c)-gammaln(length-weight-c+1))/LN2)
        factors = repetitions/2*contraction
        value = NEG
        for first in range(0, len(full), 32):
            budget.check()
            c = weight-full[first:first+32, None]-full[None, :]
            terms = counts[c]+factors[first:first+32, None]+factors[None, :]
            value = np.logaddexp2(value, np.logaddexp2.reduce(terms.ravel()))
    else:
        raise ValueError('Support method must be marginal or runs')
    budget.check()
    return float(value+max(0., -log_ratio))


class PositiveBankScreen:
    """Fixed positive regional numerators, normalized before bank minima."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, q=FIELD_ORDER,
                 distance=F(7, 20), outer_measure='invertible-toeplitz', max_seconds=60.,
                 _packet=None, _budget=None):
        self.bits = checked_q(q, construction=True)
        self.q, self.distance = q, F(distance)
        if not 0 < self.distance < F(1, 2):
            raise ValueError('Expected a distance strictly between zero and one half')
        self.z, self.alpha, self.log_ratio = float(z), float(alpha), float(log_ratio)
        self.odds = _valid_markers(marker_log_odds)
        self.outer_measure = outer_measure
        self.density = outer_log_density(q, outer_measure)
        self.forbidden = int(self.distance*OUTPUT_LENGTH)
        self.constant = -self.alpha*(self.bits+log1mexp2(-self.bits)+self.forbidden*log2(self.z))
        budget = _budget or Budget(max_seconds)
        packet = _packet or PositivePacket(self.z, self.alpha, q=q, _budget=budget)
        if (packet.z, packet.alpha, packet.q) != (self.z, self.alpha, self.q):
            raise ValueError('Cached packet parameters differ from this witness')
        self.numerators = np.array([power(packet.pair_logs(odds, _budget=budget), 64)
                                    for odds in self.odds])
        contraction_logs(self.numerators, self.log_ratio)
        budget.check()

    def regional(self, weight):
        """Every feasible full count, in ascending order, including endpoints."""
        if type(weight) is not int or not 1 <= weight <= GROUPS:
            raise ValueError('Selected weight must be in 1..8192')
        full = full_counts(PAIRS, weight)
        single, empty = weight-2*full, PAIRS-weight+full
        normalizer = LOG_FACTORIAL[PAIRS]-LOG_FACTORIAL[empty]-LOG_FACTORIAL[single]-LOG_FACTORIAL[full]
        best = np.full((len(full), 2, 2), np.inf)
        for odds, numerator in zip(self.odds, self.numerators):
            normalized = numerator[None, :, :]-normalizer[:, None, None]
            normalized -= (single*odds[0]+full*odds[1])[:, None, None]
            np.minimum(best, normalized, out=best)
        return best

    def evaluate(self, weight, *, method='marginal', optimize_potential=False,
                 max_seconds=60., max_evaluations=64, _budget=None):
        budget = _budget or Budget(max_seconds)
        budget.check()
        regions = self.regional(weight)
        ratio, evaluations = self.log_ratio, 0
        if optimize_potential:
            if not isinstance(max_evaluations, int) or not 1 <= max_evaluations <= 256:
                raise ValueError('Potential optimizer max_evaluations must lie in 1..256')
            def objective(value):
                nonlocal evaluations
                evaluations += 1
                return support_sum_logs(PAIRS, weight, regions, value, method=method, _budget=budget)
            fitted = minimize_scalar(objective, bounds=(-4096., 4096.), method='bounded',
                                     options=dict(maxiter=max_evaluations, xatol=1e-7))
            ratio = float(fitted.x)
        support = support_sum_logs(PAIRS, weight, regions, ratio, method=method, _budget=budget)
        prefix = self.alpha*self.density*weight+self.constant
        value = float(prefix+support)
        if not isfinite(value):
            raise ArithmeticError('Invalid all-profile FLOAT value')
        return dict(active_groups=weight, regional_histograms=len(regions),
                    log2_upper_screen=value, log2_support_sum_screen=support,
                    log2_outer_prefix=prefix, log2_potential_ratio=ratio,
                    potential_evaluations=evaluations, support_method=method,
                    covers_all_matching_profiles=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False)

    def box(self, ranges, *, marker_indices=(0, 0), max_seconds=60., _budget=None):
        """Convex bound for the entire clipped (r,a,b) box, not sampled points."""
        budget = _budget or Budget(max_seconds)
        ranges = validate_box(ranges)
        count = profile_count(ranges)
        if not count:
            raise ValueError('No feasible profile lies in this box')
        if len(marker_indices) != 2 or any(type(i) is not int or not 0 <= i < len(self.odds)
                                           for i in marker_indices):
            raise ValueError('Two valid fixed marker-bank indices are required')
        vertices = histogram_vertices(ranges)
        numerators = [self.numerators[i] for i in marker_indices]
        numerator = total(power(multiply(*numerators), 128)[0])
        envelopes = []
        for matching in (0, 1):
            values = []
            for r, a, b in vertices:
                budget.check()
                value = self.constant+self.alpha*self.density*float(r)
                for parity, (f, index) in enumerate(zip((a, b), marker_indices)):
                    hist = np.array((PAIRS-r+f, r-2*f, f), dtype=float)
                    mult = float((gammaln(PAIRS+1)-gammaln(hist+1).sum())/LN2)
                    value -= 128*(mult+hist[1]*self.odds[index][0]+hist[2]*self.odds[index][1])
                    if parity == matching:
                        value += mult+hist[1]
                values.append(value)
            envelopes.append(max(values))
        return dict(ranges=ranges, integer_profile_triples=count, convex_vertices=len(vertices),
                    marker_indices=marker_indices,
                    log2_upper_screen=float(numerator+min(envelopes)+log2(count)),
                    covers_entire_clipped_box=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False)

    def witness(self):
        return dict(z=str(F(self.z)), alpha=str(F(self.alpha)), field_order=self.q,
                    field_bits=self.bits, target_distance=str(self.distance),
                    forbidden_weight=self.forbidden, outer_measure=self.outer_measure,
                    marker_log_odds=self.odds, log2_potential_ratio=self.log_ratio,
                    marker_log_odds_exact=[[str(F(x)) for x in row] for row in self.odds],
                    log2_potential_ratio_exact=str(F(self.log_ratio)),
                    coefficient_method='positive-markers-log2-FLOAT', pair_swaps=True)

    def band(self, first, last, *, method='marginal', optimize_potential=False,
             max_seconds=60., max_weights=8192, threshold=-40., progress=None):
        if (type(first) is not int or type(last) is not int or not 1 <= first <= last <= GROUPS
                or not 1 <= max_weights <= GROUPS or last-first+1 > max_weights or not isfinite(threshold)):
            raise ValueError('Requested band exceeds its explicit integer weight budget')
        budget, rows, failure = Budget(max_seconds), [], None
        for weight in range(first, last+1):
            try:
                row = self.evaluate(weight, method=method, optimize_potential=optimize_potential, _budget=budget)
                rows.append(row)
                if progress:
                    progress(row)
            except (TimeoutError, ArithmeticError) as error:
                failure = str(error)
                break
        complete = len(rows) == last-first+1
        union = float(np.logaddexp2.reduce([row['log2_upper_screen'] for row in rows])) if rows else None
        return dict(requested_weight_range=[first, last], evaluated_weights=[row['active_groups'] for row in rows],
                    every_requested_weight_evaluated=complete, completed_band_geometry=complete,
                    pending_weight_range=None if complete else [first+len(rows), last],
                    failed_threshold_weights=[row['active_groups'] for row in rows if row['log2_upper_screen'] > threshold],
                    evaluated_union_log2_upper_screen=union, threshold=threshold,
                    requested_band_passes_float_threshold=complete and union <= threshold,
                    results=rows, witness=self.witness(), failure=failure,
                    elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                    covers_all_matching_profiles_at_evaluated_weights=True,
                    covers_all_even_odd_splits=True, directed=False, whole_code_certificate=False)


def terminal_tail_screen(*, q=FIELD_ORDER, distance=F(7, 20), failure_bits=40, first=None, precision=192):
    """Route-independent envelope for ALL source weights >= first.

    Fixed inner setup gives one inverse outer vector for each output. The
    target vector determines its active groups, so no support-binomial union
    occurs. The outer point cap is q^-96*(q-1)^-32 per nonzero group.
    Exact rational mark arithmetic and Arb logs avoid a floating q**k.
    This API reports a FLOAT proposal, despite using Arb for scalar stability.
    """
    bits = checked_q(q, construction=True)
    distance = F(distance)
    if (not 0 < distance < F(1, 2) or not isfinite(failure_bits) or failure_bits < 0
            or type(precision) is not int or not 128 <= precision <= 1024):
        raise ValueError('Invalid terminal target or scalar precision')
    limit = int(distance*OUTPUT_LENGTH)
    if limit < 1:
        raise ValueError('This terminal screen requires positive forbidden weight')
    mark = F(limit, (OUTPUT_LENGTH-limit)*(q-1))
    with ctx.workprec(precision):
        lg = lambda x: x.log()/arb(2).log()
        z = arb(mark.numerator)/arb(mark.denominator)
        logqm1 = lg(arb(q-1))
        ball = OUTPUT_LENGTH*lg(1+(q-1)*z)-limit*lg(z)
        loss = 96*bits+32*logqm1
        if first is None:
            first = max(1, int(float((ball-logqm1+failure_bits)/loss))+1)
            while first > 1 and ball-(first-1)*loss-logqm1 < -failure_bits:
                first -= 1
            while first <= GROUPS and not ball-first*loss-logqm1 < -failure_bits:
                first += 1
        if type(first) is not int or not 1 <= first <= GROUPS:
            raise ArithmeticError('No nonempty terminal range passes this screen target')
        value = ball-first*loss-logqm1
        numeric, enclosure = float(value), str(value)
        previous = float(ball-(first-1)*loss-logqm1) if first > 1 else None
        ball_value = float(ball)
    return dict(weight_range=[first, GROUPS], field_order=q, field_bits=bits,
                target_distance=str(distance), forbidden_weight=limit,
                mark=str(mark), log2_ball_count_screen=ball_value,
                log2_upper_screen=numeric, scalar_enclosure=enclosure,
                previous_threshold_envelope_screen=previous, requested_failure_bits=failure_bits,
                passes_float_threshold=numeric <= -failure_bits,
                covers_all_supports_in_range=True, covers_all_even_odd_splits=True,
                no_support_subset_union=True, method='inverse-image-terminal-range',
                outer_measure='invertible-toeplitz', directed=False, whole_code_certificate=False)


def source_pins():
    result = sparse_source_pins()
    directory = Path(__file__).resolve().parent
    for name in ('hd_field128_lift_tail_screen.py', 'hd_shifted_pairs_screen.py', 'hd_shifted_cover_backend.py'):
        result[name] = sha256((directory/name).read_bytes()).hexdigest()
    return result


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--terminal', action='store_true')
    parser.add_argument('--first', type=int, default=65)
    parser.add_argument('--last', type=int, default=128)
    parser.add_argument('--z', default='1/2')
    parser.add_argument('--alpha', default='1/128')
    parser.add_argument('--marker', nargs=2, type=float, action='append', metavar=('LOG_SINGLE', 'LOG_FULL'))
    parser.add_argument('--log-ratio', type=float, default=-8.)
    parser.add_argument('--method', choices=('marginal', 'runs'), default='marginal')
    parser.add_argument('--optimize-potential', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--max-weights', type=int, default=8192)
    parser.add_argument('--threshold', type=float, default=-40.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('--field-bits must be a multiple of eight in 8..128')
    if args.output and args.output.exists():
        parser.error('Output exists; choose an unused filename')
    q = 1 << args.field_bits
    if args.terminal:
        result = terminal_tail_screen(q=q, distance=F(args.distance), failure_bits=-args.threshold)
    else:
        budget = Budget(args.max_seconds)
        screen = PositiveBankScreen(float(F(args.z)), float(F(args.alpha)), args.marker or ((-4., -4.),),
                                    args.log_ratio, q=q, distance=F(args.distance), _budget=budget)
        remaining = args.max_seconds-(monotonic()-budget.started)
        budget.check()
        result = screen.band(args.first, args.last, method=args.method, optimize_potential=args.optimize_potential,
                             max_seconds=remaining, max_weights=args.max_weights, threshold=args.threshold)
    result['source_sha256'] = source_pins()
    result['pin_scope'] = 'listed source files, not all transitive imports'
    encoded = json.dumps(result, indent=2, allow_nan=False)
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            target.write(encoded+'\n')
    print(encoded)


if __name__ == '__main__':
    main()
