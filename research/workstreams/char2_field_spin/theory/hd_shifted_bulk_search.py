"""Bounded FLOAT search for disjoint unsigned global-coverage boxes.

Each witness caches its positive regional matrix products once. Box checks
then use only affine forms, cached log-factorials, and clipped vertices.
Accepted boxes are explicit proof proposals, not directed certificates.
The returned partition also retains every boundary and unfinished box.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import heapq
from math import lgamma, log, log2, prod

import numpy as np

from hd_grouped_packet_screen import LIMIT
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_unoriented_screen import positive_pair_matrix
from hd_shifted_cover_backend import CachedBound, histogram_vertices, profile_count


PAIRS = 4096
LOG_FACTORIAL = np.array([lgamma(n+1)/log(2) for n in range(PAIRS+1)])
LOG_DENSITY = 2048-128*log2(255)

# A targeted correction to the initial bank's fixed-alpha family. These
# overlapping boxes are diagnostic witnesses, not a disjoint cover plan.
ALPHA_MARK = .011103584869424775
ALPHA_EXPONENT = .5357886737301887
ALPHA_CHECKPOINTS = (
    (((5120, 5247), (1024, 1278), (1024, 1231)),
     (7.905242031272326, 9.62627537577853), -244.6021153635386),
    (((5120, 5151), (1024, 1055), (1024, 1055)),
     (11.095348221791237, 12.764697527220385), -5440.182220670806),
    (((5120, 5247), (1152, 1407), (1152, 1407)),
     (6.907696045558403, 8.905543123712771), -31487.133480863296),
    (((5248, 5375), (1280, 1535), (1280, 1535)),
     (6.780383317529136, 8.905285738995072), -67094.99943455844),
    (((5120, 5375), (1280, 1535), (1280, 1535)),
     (6.283482967970975, 8.617029057088011), -52742.0022156348),
    (((5120, 5247), (1024, 1087), (1024, 1087)),
     (10.086809951318948, 11.760985326348134), -4843.336308509502),
)


def volume(box):
    return prod(hi-lo+1 for lo, hi in box)


def split_box(box):
    widths = tuple(hi-lo for lo, hi in box)
    axis = max(range(3), key=widths.__getitem__)
    if widths[axis] == 0:
        raise ValueError('Cannot split a singleton box')
    lo, hi = box[axis]
    middle = (lo+hi)//2
    left, right = list(box), list(box)
    left[axis] = lo, middle
    right[axis] = middle+1, hi
    return tuple(left), tuple(right)


@lru_cache(None)
def box_geometry(box):
    """Return clipped vertices and both one-matching entropy terms."""
    vertices = histogram_vertices(box)
    if not vertices or not profile_count(box):
        raise ValueError('No feasible profile in box')
    corners = np.array(vertices, dtype=float)
    r, a, b = corners.T
    hist0 = np.stack((PAIRS-r+a, r-2*a, a), axis=1)
    hist1 = np.stack((PAIRS-r+b, r-2*b, b), axis=1)
    factorial = np.vectorize(lambda n: lgamma(n+1)/log(2))
    first = LOG_FACTORIAL[PAIRS]-factorial(hist0).sum(axis=1)
    second = LOG_FACTORIAL[PAIRS]-factorial(hist1).sum(axis=1)
    return corners, np.stack((-127*first-128*second, -128*first-127*second))


@dataclass(frozen=True)
class PositiveWitness:
    z: float
    alpha: float
    odds: tuple
    centre: tuple
    intercept: float
    slopes: tuple
    log2_numerator: float

    @classmethod
    def create(cls, z, alpha, odds, centre=(0, 0, 0)):
        odds = tuple(tuple(map(float, row)) for row in odds)
        matrices = {}
        regions = []
        for row in odds:
            if row not in matrices:
                positive = positive_pair_matrix(z, alpha, row)
                if not np.all(positive > 0) or not np.all(np.isfinite(positive)):
                    raise ArithmeticError('Invalid positive matrix evaluation')
                matrices[row] = power(tuple(map(tuple, np.log2(positive))), 64)
            regions.append(matrices[row])
        value = power(multiply(regions[0], regions[1]), 128)
        numerator = float(total(value[0]))
        intercept = numerator-alpha*(log2(255)+LIMIT*log2(z))
        (s0, f0), (s1, f1) = odds
        slopes = (alpha*LOG_DENSITY+1-128*(s0+s1),
                  -2+256*s0-128*f0, 256*s1-128*f1)
        return cls(float(z), float(alpha), odds, tuple(centre), intercept, slopes, numerator)

    def evaluate(self, box):
        corners, entropy = box_geometry(box)
        slopes = np.array((self.slopes, np.asarray(self.slopes)+(0, 2, -2)))
        values = self.intercept+slopes@corners.T+entropy
        matching = int(np.argmin(np.max(values, axis=1)))
        index = int(np.argmax(values[matching]))
        maximum = float(values[matching, index])
        return maximum+log2(profile_count(box)), tuple(corners[index])

    def serializable(self):
        return dict(mark=self.z, alpha=self.alpha, marker_log_odds=[list(row) for row in self.odds],
                    proposal_centre=list(self.centre), log2_positive_numerator=self.log2_numerator,
                    cached_intercept=self.intercept, cached_slopes=list(self.slopes))


def seeded_witness(box, tier_center):
    """A density-scaled proposal; acceptance always checks the whole box."""
    centre = tuple(np.mean(box_geometry(box)[0], axis=0))
    r, a, b = centre
    z = .23149023644229288**(tier_center/2048)
    alpha = .11372673628845778
    odds = []
    for full in (a, b):
        empty, singleton = PAIRS-r+full, r-2*full
        odds.append((log2(max(1, singleton)/max(1, empty))+.5182296355378473,
                     log2(max(1, full)/max(1, empty))+.9361653072364948))
    return PositiveWitness.create(z, alpha, odds, centre)


class AdaptiveAlphaProposer:
    """Optional bounded witness refinement; it never certifies a box itself.

    First refine four regional marker odds at a fixed (z,alpha). Try the
    validated alpha and then a finite alpha family. Only a remaining
    failure permits a short weight-marker refinement. All trial witnesses
    use the same whole-box envelope as the ordinary search.
    """
    witness_type = PositiveWitness

    def __init__(self, *, max_evaluations=12000, marker_maxfev=200,
                 weight_maxfev=32, alpha_family=(.25, .5, .75, 1.)):
        if (max_evaluations < 1 or marker_maxfev < 1 or weight_maxfev < 0
                or not alpha_family or any(not 0 < a <= 1 for a in alpha_family)):
            raise ValueError('Positive work limits and admissible alpha candidates required')
        self.max_evaluations = int(max_evaluations)
        self.marker_maxfev = int(marker_maxfev)
        self.weight_maxfev = int(weight_maxfev)
        self.alpha_family = tuple(dict.fromkeys((ALPHA_EXPONENT, *map(float, alpha_family))))
        self.evaluations = self.marker_optimizations = self.weight_optimizations = 0
        self.box_attempts = 0

    @property
    def exhausted(self):
        return self.evaluations >= self.max_evaluations

    def create(self, z, alpha, odds, centre=(0, 0, 0)):
        if self.exhausted:
            return None
        self.evaluations += 1
        return self.witness_type.create(z, alpha, odds, centre)

    def initial_bank(self):
        result = []
        for _, odds, _ in ALPHA_CHECKPOINTS:
            if self.exhausted:
                break
            result.append(self.create(ALPHA_MARK, ALPHA_EXPONENT, (odds, odds)))
        return result

    @staticmethod
    def marker_seed(box, alpha):
        centre = np.mean(box_geometry(box)[0], axis=0)
        r, a, b = centre
        # Translate the validated interior witness by histogram log odds.
        offset = np.array((6.22715071206585-log2(10),
                           8.584298238117446-log2(5)))*(alpha/ALPHA_EXPONENT)
        rows = []
        for full in (a, b):
            empty, singleton = PAIRS-r+full, r-2*full
            rows.append(np.log2(np.maximum((singleton, full), 1)/max(empty, 1))+offset)
        return np.clip(np.array(rows).ravel(), -15, 13.5), tuple(centre)

    def propose(self, box, tier, threshold):
        from scipy.optimize import minimize

        self.box_attempts += 1
        best, score = None, float('inf')

        class SufficientWitness(Exception):
            pass

        class WorkLimit(Exception):
            pass

        def objective(parameters, alpha, centre, vary_weight=False):
            nonlocal best, score
            if self.exhausted:
                raise WorkLimit
            # Every attempted numerical objective evaluation consumes budget,
            # including rejected marker proposals.
            self.evaluations += 1
            x = np.asarray(parameters)
            log_z, odds = (x[0], x[1:]) if vary_weight else (log2(ALPHA_MARK), x)
            if (not np.all(np.isfinite(x)) or not -14 <= log_z <= -.01
                    or np.min(odds) < -16 or np.max(odds) > 14):
                return 1e15
            try:
                mark = 2**log_z if vary_weight else ALPHA_MARK
                witness = self.witness_type.create(mark, alpha, odds.reshape(2, 2), centre)
                value = witness.evaluate(box)[0]
            except (ArithmeticError, FloatingPointError, OverflowError):
                return 1e15
            if value < score:
                best, score = witness, value
            if score <= threshold:
                raise SufficientWitness
            return value

        try:
            for alpha in self.alpha_family:
                if self.exhausted:
                    break
                seed, centre = self.marker_seed(box, alpha)
                self.marker_optimizations += 1
                minimize(objective, seed, args=(alpha, centre), method='Nelder-Mead',
                         options=dict(maxfev=min(self.marker_maxfev,
                                               self.max_evaluations-self.evaluations),
                                      xatol=1e-6, fatol=1e-4))
            if best is not None and self.weight_maxfev and not self.exhausted:
                self.weight_optimizations += 1
                seed = np.r_[log2(best.z), np.array(best.odds).ravel()]
                minimize(objective, seed, args=(best.alpha, best.centre, True), method='Nelder-Mead',
                         options=dict(maxfev=min(self.weight_maxfev,
                                               self.max_evaluations-self.evaluations),
                                      xatol=1e-5, fatol=1e-3))
        except (SufficientWitness, WorkLimit):
            pass
        return best

    def statistics(self):
        return dict(max_evaluations=self.max_evaluations, evaluations=self.evaluations,
                    marker_maxfev=self.marker_maxfev, weight_maxfev=self.weight_maxfev,
                    alpha_family=list(self.alpha_family), box_attempts=self.box_attempts,
                    marker_optimizations=self.marker_optimizations,
                    weight_optimizations=self.weight_optimizations)


def chunks(low, high, width):
    return tuple((first, min(high, first+width-1)) for first in range(low, high+1, width))


def initial_partition(first=4096, last=7648, *, slab_width=256, full_width=256, pairs=PAIRS):
    """Disjoint interiors plus four boundary rectangles per weight slab.

    Their union contains every histogram-admissible triple in the requested
    weight range. Only interior boxes have nonnegative histograms at every
    corner. The search handles boundary rectangles using clipped vertices.
    """
    if not 1 <= first <= last < 2*pairs or min(slab_width, full_width) < 1:
        raise ValueError('Invalid range or partition widths')
    interiors, boundaries = [], []
    for low, high in chunks(first, last, slab_width):
        outer_low, outer_high = max(0, low-pairs), high//2
        inner_low, inner_high = max(0, high-pairs), low//2
        if inner_low > inner_high:
            boundaries.append(((low, high), (outer_low, outer_high), (outer_low, outer_high)))
            continue
        for a in chunks(inner_low, inner_high, full_width):
            for b in chunks(inner_low, inner_high, full_width):
                interiors.append((((low, high), a, b), (low+high)/2))
        outer = (outer_low, outer_high)
        inner = (inner_low, inner_high)
        edges = ((outer_low, inner_low-1), (inner_high+1, outer_high))
        for edge in edges:
            if edge[0] <= edge[1]:
                boundaries.append(((low, high), edge, outer))
                boundaries.append(((low, high), inner, edge))
    return interiors, boundaries


def best_witness(bank, box):
    if not bank:
        return None, float('inf'), None
    corners, entropy = box_geometry(box)
    intercepts = np.array([w.intercept for w in bank])
    base = np.array([w.slopes for w in bank])
    slopes = np.stack((base, base+(0, 2, -2)), axis=1)
    values = intercepts[:, None, None]+slopes@corners.T+entropy[None, :, :]
    maxima = np.min(np.max(values, axis=2), axis=1)+log2(profile_count(box))
    index = int(np.argmin(maxima))
    matching = int(np.argmin(np.max(values[index], axis=1)))
    worst = tuple(corners[int(np.argmax(values[index, matching]))])
    return index, float(maxima[index]), worst


def search(*, first=4096, last=7648, proposals=400, threshold=-96.,
           slab_width=256, full_width=256, minimum_width=16, progress=None, resume=None,
           adaptive=None):
    if proposals < 1 or minimum_width < 1:
        raise ValueError('Positive finite proposal budget and minimum width required')
    bank, accepted = [], []
    if resume is None:
        initial, boundaries = initial_partition(first, last, slab_width=slab_width, full_width=full_width)
    else:
        source_first, source_last = resume['requested_weight_range']
        if not source_first <= first <= last <= source_last:
            raise ValueError('Continuation must stay inside the saved domain')
        bank = [PositiveWitness.create(row['mark'], row['alpha'], row['marker_log_odds'],
                                       row['proposal_centre']) for row in resume['witnesses']]
        initial, boundaries = [], []
        source_width = resume['slab_width']
        for inherited, rows in ((True, resume['accepted_boxes']), (False, resume['unresolved_boxes'])):
            for row in rows:
                (lo, hi), a, b = row['ranges']
                box = ((max(first, lo), min(last, hi)), tuple(a), tuple(b))
                if box[0][0] > box[0][1] or not profile_count(box):
                    continue
                if inherited:
                    score, worst = bank[row['witness_id']].evaluate(box)
                    accepted.append(dict(ranges=box, witness_id=row['witness_id'], inherited=True,
                                         union_log2_upper_screen=score, worst_corner=worst,
                                         integer_triples=profile_count(box)))
                else:
                    slab = source_first+((box[0][0]-source_first)//source_width)*source_width
                    tier = (slab+min(source_last, slab+source_width-1))/2
                    initial.append((box, tier))
    if adaptive is not None:
        bank.extend(adaptive.initial_bank())
    pending, serial = [], 0
    def enqueue(box, tier, depth=0):
        nonlocal serial
        if not profile_count(box):
            return
        heapq.heappush(pending, (depth, -box[0][1], -volume(box), serial, box, tier))
        serial += 1
    for box, tier in initial:
        enqueue(box, tier)
    for box in boundaries:
        enqueue(box, sum(box[0])/2)
    unresolved = []
    used, reported = 0, 0
    while pending:
        depth, _, _, _, box, tier = heapq.heappop(pending)
        witness_id, score, worst = best_witness(bank, box)
        if score > threshold:
            if used >= proposals or adaptive is not None and adaptive.exhausted:
                reason = 'proposal-budget' if used >= proposals else 'adaptive-evaluation-budget'
                unresolved.append(dict(ranges=box, reason=reason, best_log2_upper_screen=score))
                continue
            candidate = (seeded_witness(box, tier) if adaptive is None
                         else adaptive.propose(box, tier, threshold))
            if candidate is not None:
                bank.append(candidate)
            used += 1
            witness_id, score, worst = best_witness(bank, box)
        if score <= threshold:
            accepted.append(dict(ranges=box, witness_id=witness_id,
                                 union_log2_upper_screen=score, worst_corner=worst,
                                 integer_triples=profile_count(box), inherited=False))
        elif max(hi-lo+1 for lo, hi in box) <= minimum_width:
            unresolved.append(dict(ranges=box, reason='minimum-width', best_log2_upper_screen=score))
        else:
            for child in split_box(box):
                enqueue(child, tier, depth+1)
        if progress and used//40 > reported:
            reported = used//40
            progress(dict(stage='bulk-search', proposals=used, accepted_boxes=len(accepted),
                          accepted_triples=sum(row['integer_triples'] for row in accepted),
                          pending_boxes=len(pending)))
    accepted_count = sum(row['integer_triples'] for row in accepted)
    domain = profile_count(((first, last), (0, PAIRS), (0, PAIRS)))
    maximum = max((row['union_log2_upper_screen'] for row in accepted), default=-float('inf'))
    union = (maximum+log2(sum(2**(row['union_log2_upper_screen']-maximum) for row in accepted))
             if accepted else None)
    result = dict(profile_coordinates=['r', 'a', 'b'], all_even_odd_splits=True,
                requested_weight_range=[first, last], proposals_used=used, proposal_budget=proposals,
                threshold=threshold, slab_width=slab_width, full_width=full_width,
                minimum_width=minimum_width, accepted_boxes=accepted, unresolved_boxes=unresolved,
                inherited_accepted_boxes=sum(row.get('inherited', False) for row in accepted),
                witnesses=[w.serializable() for w in bank],
                accepted_integer_triples=accepted_count,
                feasible_domain_triples=domain,
                accepted_domain_fraction=accepted_count/domain,
                union_log2_upper_screen=union, directed=False, whole_code_certificate=False)
    if adaptive is not None:
        result['adaptive_proposer'] = adaptive.statistics()
    return result


def compact_plan(result):
    """Retain explicit leaves and only their used numerical witnesses."""
    used = sorted({row['witness_id'] for row in result['accepted_boxes']})
    mapping = {old: new for new, old in enumerate(used)}
    fields = ('mark', 'alpha', 'marker_log_odds', 'proposal_centre')
    witnesses = [{key: result['witnesses'][old][key] for key in fields} for old in used]
    accepted = []
    for row in result['accepted_boxes']:
        new = dict(row, witness_id=mapping[row['witness_id']])
        new.pop('worst_corner', None)
        accepted.append(new)
    return dict(result, witnesses=witnesses, accepted_boxes=accepted)


def alpha_checkpoint_bounds():
    """Replay six fixed proposals without running a witness search."""
    result = []
    for box, odds, _ in ALPHA_CHECKPOINTS:
        witness = PositiveWitness.create(ALPHA_MARK, ALPHA_EXPONENT, (odds, odds))
        bound, worst = witness.evaluate(box)
        result.append(dict(ranges=box, integer_triples=profile_count(box),
                           union_log2_upper_screen=bound, worst_vertex=worst,
                           witness=witness.serializable(), directed=False,
                           whole_code_certificate=False))
    return result


def directed_alpha_checkpoints(precision=256):
    """Replay the three selected correction boxes with exact input markers."""
    result = []
    for index in (0, 4, 2):
        box, odds, _ = ALPHA_CHECKPOINTS[index]
        row = tuple(Fraction(2**value).limit_denominator(10**9) for value in odds)
        witness = CachedBound(Fraction(ALPHA_MARK), Fraction(ALPHA_EXPONENT),
                              (row, row), precision=precision)
        result.append(dict(checkpoint_index=index, witness=witness.witness(), **witness.receipt(box)))
    return result


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, default=4096)
    parser.add_argument('--last', type=int, default=7648)
    parser.add_argument('--proposals', type=int, default=400)
    parser.add_argument('--slab-width', type=int, default=256)
    parser.add_argument('--full-width', type=int, default=256)
    parser.add_argument('--compact', action='store_true', help='Omit unused witnesses and cached arithmetic')
    parser.add_argument('--resume', type=str, help='Read a saved explicit partition and refine its unresolved leaves')
    parser.add_argument('--alpha-checkpoints', action='store_true', help='Replay six fixed alpha-correction boxes only')
    parser.add_argument('--precision', type=int, help='With --alpha-checkpoints, replay three selected boxes using Arb')
    parser.add_argument('--adaptive-alpha', action='store_true', help='Use the optional bounded alpha/marker proposer')
    parser.add_argument('--max-evaluations', type=int, default=12000, help='Adaptive objective-call budget')
    parser.add_argument('--marker-maxfev', type=int, default=200, help='Per-alpha marker objective-call limit')
    parser.add_argument('--weight-maxfev', type=int, default=32, help='Weight-marker refinement call limit; zero disables')
    args = parser.parse_args()
    if args.precision is not None and not args.alpha_checkpoints:
        parser.error('--precision is only for the fixed --alpha-checkpoints replay')
    if args.alpha_checkpoints:
        result = (alpha_checkpoint_bounds() if args.precision is None
                  else directed_alpha_checkpoints(args.precision))
        print(json.dumps(dict(stage='alpha-checkpoints', results=result), allow_nan=False))
        raise SystemExit(0)
    if args.resume:
        with open(args.resume, encoding='utf-8') as source:
            saved = json.load(source)
    else:
        saved = None
    adaptive = (AdaptiveAlphaProposer(max_evaluations=args.max_evaluations,
                                     marker_maxfev=args.marker_maxfev,
                                     weight_maxfev=args.weight_maxfev) if args.adaptive_alpha else None)
    result = search(first=args.first, last=args.last, proposals=args.proposals, resume=saved, adaptive=adaptive,
                    slab_width=args.slab_width, full_width=args.full_width,
                    progress=lambda item: print(json.dumps(item), flush=True))
    if args.compact:
        result = compact_plan(result)
    print(json.dumps(dict(stage='bulk-result', result=result), allow_nan=False), flush=True)
