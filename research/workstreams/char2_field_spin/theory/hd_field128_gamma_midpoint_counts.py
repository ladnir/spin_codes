"""Joint matching-count caps using the proved midpoint maximum.

Two one-dimensional exact binomial arrays replace a separate flow sum for
every profile. Counts are for indexed words. FLOAT logarithms are for
proposal screens only; count() uses integers and checks exact division.
"""

from hashlib import sha256
from math import comb, log2
from pathlib import Path

import numpy as np

import hd_field128_gamma_exact_profile_counts as exact
import hd_field128_gamma_near_face_counts as near


NOTE = Path(__file__).resolve().parent.parent/'HD_FIELD128_MIDPOINT_COUNTS.md'


def source_pins():
    pins = {**near.source_pins(), **exact.source_pins()}
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    pins['../'+NOTE.name] = sha256(NOTE.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


class MidpointCounts:
    """All-profile upper counts, with O(r) exact preparation and storage."""

    def __init__(self, pairs, weight, *, max_singletons=64, budget=None):
        budget = budget or near.provider.old.Budget(60.)
        self.base = near.NearFaceCounts(pairs, weight, max_singletons=max_singletons, budget=budget)
        self.pairs, self.weight, self.full = pairs, weight, self.base.full
        self.shift = max(0, weight-pairs)
        self.r = min(weight, 2*pairs-weight)
        self.low, self.high = self.r//2, (self.r+1)//2
        self.denominator = (pairs-self.low)*self.low*(pairs-self.high)*self.high
        left, right = [], []
        for a in range(self.low+1):
            if not a % 64:
                budget.check()
            left.append(comb(pairs-self.low, self.high-a)*comb(self.low, a))
            right.append(comb(pairs-self.high, self.low-a)*comb(self.high, a))
        self.left, self.right = tuple(left), tuple(right)
        self._left_logs = np.array([log2(v) for v in left])
        self._right_logs = np.array([log2(v) for v in right])
        self._constant = log2(pairs)-log2(self.denominator) if self.r >= 2 else 0.
        budget.check()

    def midpoint_upper(self, first, second):
        if type(first) is not int or type(second) is not int:
            raise ValueError('Integer full-pair counts required')
        if first not in self.full or second not in self.full:
            return 0
        a, b = max(first, second)-self.shift, min(first, second)-self.shift
        if self.r <= 1:
            return 1 if self.r == 0 else 2*self.pairs
        if self.r == self.pairs and a == b == 0:
            return 2
        tau = ((self.low-a)*(self.pairs-self.r+b)*self.high
               + a*(self.high-b)*(self.pairs-self.high))
        value, remainder = divmod(self.pairs*tau*self.left[a]*self.right[b], self.denominator)
        if remainder:
            raise ArithmeticError('Nonintegral midpoint flow count')
        return (self.r-2*a+1)*value

    def count(self, first, second):
        return min(self.base.count(first, second), self.midpoint_upper(first, second))

    def log_counts(self, first, second, *, budget=None):
        baseline = self.base.log_counts(first, second, budget=budget)
        first, second = np.broadcast_arrays(np.asarray(first), np.asarray(second))
        a, b = np.maximum(first, second)-self.shift, np.minimum(first, second)-self.shift
        if self.r <= 1:
            return np.minimum(baseline, 0. if self.r == 0 else log2(2*self.pairs))
        tau = ((self.low-a)*(self.pairs-self.r+b)*self.high
               + a*(self.high-b)*(self.pairs-self.high))
        with np.errstate(divide='ignore'):
            values = (self._constant+np.log2(tau)+np.log2(self.r-2*a+1)
                      + self._left_logs[a]+self._right_logs[b])
        if self.r == self.pairs:
            values = np.where((a == 0) & (b == 0), 1., values)
        return np.minimum(baseline, values)

    def diagnostics(self):
        return dict(method='min(previous valid profile cap, (s+1) times exact midpoint flow count)',
                    pairs=self.pairs, active_groups=self.weight, canonical_weight=self.r,
                    maximum_refined_singletons=self.base.max_singletons,
                    all_profiles_retained=True, integer_count_path=True,
                    maximum_factor_over_exact_joint_count=self.r+1,
                    alternating_exception='r=P,a=b=0 has exact count2',
                    prepared_integer_array_entries=len(self.left)+len(self.right),
                    no_profile_pair_table=True, source_proof_note='../'+NOTE.name)


def support_sum(pairs, weight, regions, *, repetitions=256, chunk=32,
                counts=None, budget=None, max_singletons=64):
    """Complete FLOAT transfer, for proposing proof parameters only."""
    transfer = near.provider.profiles.transfer
    full = transfer.tail.full_counts(pairs, weight)
    regions = np.asarray(regions)
    if (regions.shape != (len(full), 2, 2) or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or type(chunk) is not int or not 1 <= chunk <= 128):
        raise ValueError('Invalid regional matrices, repetitions, or chunk size')
    budget = budget or near.provider.old.Budget(60.)
    counts = counts or MidpointCounts(pairs, weight, max_singletons=max_singletons, budget=budget)
    if type(counts) is not MidpointCounts or (counts.pairs, counts.weight) != (pairs, weight):
        raise ValueError('Expected matching midpoint count object')
    total, largest = -np.inf, None
    for start in range(0, len(full), chunk):
        budget.check()
        logs = counts.log_counts(full[start:start+chunk, None], full[None, :], budget=budget)
        product = transfer.multiply(regions[start:start+chunk, None], regions[None])
        repeated = transfer.power(product, repetitions//2)
        output = np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
        terms = logs+output
        position = np.unravel_index(np.argmax(terms), terms.shape)
        candidate = float(terms[position])
        if largest is None or candidate > largest['log2_term']:
            even, odd = int(full[start+position[0]]), int(full[position[1]])
            largest = dict(even_full=even, odd_full=odd, runs=weight-even-odd,
                           log2_term=candidate, log2_count_upper=float(logs[position]),
                           log2_transfer=float(output[position]))
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(terms.ravel())))
    budget.check()
    return dict(log2_support_sum_screen=total, largest_profile=largest,
                profile_pairs_evaluated=len(full)**2, common_potential_used=False,
                count_provenance=counts.diagnostics(), all_profiles_retained=True,
                arithmetic='FLOAT transfer; proved midpoint integer caps',
                directed=False, whole_code_certificate=False)
