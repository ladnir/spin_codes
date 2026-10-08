"""Floating complete interval cover of active-row counts, not a certificate.

Marks y,tau remain fixed inside each interval. Positive regional coefficients
give R_A <= R_b entrywise for A<=b. The source factor is maximized explicitly
over every integer in the interval. All interval unions are then added.
"""

import argparse
from dataclasses import dataclass
from functools import lru_cache
import json
import math

import numpy as np
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose
from reduced_distance_certificate import log_polynomial_product
from superregular_return_screen import local_matrix


# Values returned by the bounded two-start search; logs are natural logs.
# Their y values must be evaluated at their anchor, not at a bin endpoint.
ANCHORS = {
    1: (7.268964250212447, 5.541263545158426),
    2: (7.265526562087674, 5.541263545158426),
    4: (7.263803401017213, 5.541263545158426),
    8: (7.419376729468748, 5.270829056060071),
    16: (7.183071336951522, 2.947651636355612),
    32: (7.120435500608540, 2.258220251564710),
    64: (7.062633053108466, 1.991286947926587),
    128: (7.024426560335106, 1.907863435406201),
    256: (7.028066385865181, 1.923294185949148),
    512: (7.0745259685264505, 2.0640527995049034),
    600: (7.087584031696039, 2.125052179062638),
}


# Compact replay fixture: (lower, upper, fixed-mark anchor). These are not
# probability estimates or certificates; replay evaluates every interval.
FROZEN_INTERVALS = (
    (1,1,1), (2,2,2), (3,3,3), (4,4,4), (5,6,6), (7,9,11), (10,13,16),
    (14,16,16), (17,19,16), (20,23,23), (24,27,23), (28,32,32), (33,37,32),
    (38,41,45), (42,45,45), (46,49,45), (50,53,45), (54,59,64), (60,65,64),
    (66,71,64), (72,76,64), (77,80,91), (81,84,91), (85,88,91), (89,92,91),
    (93,96,91), (97,100,91), (101,104,91), (105,107,91), (108,113,128),
    (114,119,128), (120,125,128), (126,130,128), (131,136,128), (137,141,128),
    (142,144,128), (145,147,128), (148,150,128), (151,152,128), (153,156,181),
    (157,160,181), (161,168,181), (169,176,181), (177,184,181), (185,192,181),
    (193,200,181), (201,204,181), (205,208,181), (209,212,181), (213,215,181),
    (216,221,256), (222,227,256), (228,233,256), (234,238,256), (239,249,256),
    (250,260,256), (261,271,256), (272,282,256), (283,288,256), (289,293,256),
    (294,299,256), (300,304,256), (305,312,362), (313,320,362), (321,336,362),
    (337,352,362), (353,367,362), (368,383,362), (384,399,362), (400,415,362),
    (416,423,362), (424,430,362), (431,443,512), (444,456,512), (457,481,512),
    (482,507,512), (508,532,512), (533,554,554), (555,576,554), (577,600,600),
)


def marks(anchor, point, output_length):
    if anchor <= 0 or output_length <= 0:
        raise ValueError('positive anchor and output length required')
    logtilt, logtau = point
    logy = -math.exp(logtilt)*anchor/output_length
    return logy, float(logtau)


def interpolated_anchors(anchors):
    """One extra fixed witness per geometric gap; no optimization or claim."""
    result = dict(anchors)
    ordered = sorted(anchors)
    for left, right in zip(ordered, ordered[1:]):
        middle = round(math.sqrt(left*right))
        if not left < middle < right:
            continue
        fraction = math.log(middle/left)/math.log(right/left)
        result[middle] = tuple(a+(b-a)*fraction for a, b in zip(anchors[left], anchors[right]))
    return result


def log_coefficients(y, tau, rows, degree, state=8):
    """The existing positive log polynomial, computed once for fixed marks."""
    if rows <= 0 or rows % 16 or not 0 <= degree <= rows:
        raise ValueError('invalid coefficient geometry')
    local = np.asarray([local_matrix(w, y, tau, state) for w in range(min(16, degree)+1)], dtype=float)
    with np.errstate(divide='ignore'):
        coefficients = np.log(local)
    for weight in range(len(local)):
        coefficients[weight] += math.log(math.comb(16, weight))
    result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
    power = rows//16
    while power:
        if power & 1:
            result = log_polynomial_product(result, coefficients, degree)
        power >>= 1
        if power:
            coefficients = log_polynomial_product(coefficients, coefficients, degree)
    return result


def regional_from_coefficients(coefficients, active, rows, logtau):
    if not 0 <= active <= rows or active >= len(coefficients):
        raise ValueError('active count exceeds available coefficients')
    weights = np.asarray([log_choose(active, h)-log_choose(rows, h)+h*logtau
                          for h in range(active+1)])
    return logsumexp(coefficients[:active+1]+weights[:, None, None], axis=0)


def log_terminal(matrix, length):
    values = np.asarray((0., -math.inf))
    for _ in range(length):
        values = logsumexp(values[:, None]+matrix, axis=0)
    return float(logsumexp(values))


@dataclass(frozen=True)
class Interval:
    lower: int
    upper: int
    anchor: int
    per_class_margin_bits: float
    union_margin_bits: float
    source_maximizer: int


class FixedMarks:
    def __init__(self, anchor, point, degree, k_log=20, length=256, dimension=176, state=8):
        self.anchor = anchor
        self.g = geometry(k_log, length, dimension)
        self.logy, self.logtau = marks(anchor, point, self.g['n'])
        self.coefficients = log_coefficients(math.exp(self.logy), math.exp(self.logtau),
                                             self.g['rows'], degree, state)
        parity = length-dimension
        self.source = tuple(log_choose(self.g['rows'], a)-(parity*a+1)*self.logtau
                            for a in range(degree+1))

    @lru_cache(maxsize=None)
    def terminal(self, upper):
        matrix = regional_from_coefficients(self.coefficients, upper, self.g['rows'], self.logtau)
        return log_terminal(matrix, self.g['length'])

    def interval(self, lower, upper):
        if not 1 <= lower <= upper < len(self.source):
            raise ValueError('invalid positive active-row interval')
        maximizing = max(range(lower, upper+1), key=self.source.__getitem__)
        log_bound = self.source[maximizing]-self.g['threshold']*self.logy+self.terminal(upper)
        margin = -log_bound/math.log(2)
        return Interval(lower, upper, self.anchor, margin,
                        margin-math.log2(upper-lower+1), maximizing)


def anchor_zones(lower, upper, anchors):
    """Partition integers by nearest log-scale anchor; not a proof step."""
    if not 1 <= lower <= upper or not anchors:
        raise ValueError('invalid range or empty anchor collection')
    zones = []
    start = lower
    previous = None
    for active in range(lower, upper+1):
        selected = min(anchors, key=lambda anchor: (abs(math.log(active/anchor)), anchor))
        if previous is not None and selected != previous:
            zones.append((start, active-1, previous))
            start = active
        previous = selected
    zones.append((start, upper, previous))
    return zones


def subdivide(fixed, lower, upper, minimum_margin):
    result = fixed.interval(lower, upper)
    if result.union_margin_bits >= minimum_margin:
        return [result]
    if lower == upper:
        raise ArithmeticError(f'active={lower} has insufficient floating margin: {result.union_margin_bits}')
    middle = (lower+upper)//2
    return (subdivide(fixed, lower, middle, minimum_margin)
            + subdivide(fixed, middle+1, upper, minimum_margin))


def subdivide_candidates(candidates, lower, upper, minimum_margin):
    """Every candidate is independently valid; choose the strongest one."""
    result = max((fixed.interval(lower, upper) for fixed in candidates),
                 key=lambda item: item.union_margin_bits)
    if result.union_margin_bits >= minimum_margin:
        return [result]
    if lower == upper:
        raise ArithmeticError(f'active={lower} has insufficient floating margin: {result.union_margin_bits}')
    middle = (lower+upper)//2
    return (subdivide_candidates(candidates, lower, middle, minimum_margin)
            + subdivide_candidates(candidates, middle+1, upper, minimum_margin))


def cover(lower=1, upper=600, minimum_margin=100., anchors=None, progress=None):
    anchors = interpolated_anchors(ANCHORS) if anchors is None else anchors
    ordered = sorted(anchors)
    intervals = []
    for lo, hi, anchor in anchor_zones(lower, upper, anchors):
        index = ordered.index(anchor)
        neighboring = ordered[max(0, index-1):min(len(ordered), index+2)]
        candidates = [FixedMarks(a, anchors[a], hi) for a in neighboring]
        interval_part = subdivide_candidates(candidates, lo, hi, minimum_margin)
        intervals.extend(interval_part)
        if progress:
            progress(dict(zone=[lo, hi], anchor=anchor, intervals=len(interval_part),
                          minimum_margin_bits=min(x.union_margin_bits for x in interval_part)))
    return summarize(intervals, lower, upper)


def summarize(intervals, lower, upper):
    if intervals[0].lower != lower or intervals[-1].upper != upper or any(
            left.upper+1 != right.lower for left, right in zip(intervals, intervals[1:])):
        raise ArithmeticError('interval cover has a gap or overlap')
    union_margin = -float(logsumexp([-item.union_margin_bits*math.log(2) for item in intervals]))/math.log(2)
    return dict(range=[lower, upper], interval_count=len(intervals),
                minimum_interval_margin_bits=min(x.union_margin_bits for x in intervals),
                total_union_margin_bits=union_margin, rigorous=False, whole_code=False,
                intervals=[item.__dict__ for item in intervals])


def replay_frozen():
    anchors = interpolated_anchors(ANCHORS)
    endpoint = {anchor: max(hi for _, hi, a in FROZEN_INTERVALS if a == anchor)
                for _, _, anchor in FROZEN_INTERVALS}
    fixed = {anchor: FixedMarks(anchor, anchors[anchor], hi) for anchor, hi in endpoint.items()}
    intervals = [fixed[anchor].interval(lo, hi) for lo, hi, anchor in FROZEN_INTERVALS]
    return summarize(intervals, 1, 600)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lower', type=int, default=1)
    parser.add_argument('--upper', type=int, default=600)
    parser.add_argument('--minimum-margin', type=float, default=100.)
    parser.add_argument('--show-intervals', action='store_true')
    parser.add_argument('--frozen', action='store_true', help='replay the literal complete 1..600 fixture')
    args = parser.parse_args()
    result = replay_frozen() if args.frozen else cover(
        args.lower, args.upper, args.minimum_margin,
        progress=lambda record: print(json.dumps(record), flush=True))
    if not args.show_intervals:
        del result['intervals']
    print(json.dumps(result), flush=True)
