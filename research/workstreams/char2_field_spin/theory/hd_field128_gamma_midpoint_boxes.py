"""Exact rectangle sums of the proved midpoint support-count caps.

The cap is not the exact joint count. This module sums the raw midpoint
cap, not the pointwise minimum with another cap. Six one-dimensional
integer prefixes suffice; no matching-profile pair table is built.
"""

from hashlib import sha256
from pathlib import Path

import hd_field128_gamma_midpoint_counts as midpoint


NOTE = Path(__file__).resolve().parent.parent/'HD_FIELD128_MIDPOINT_BOXES.md'


def source_pins():
    pins = midpoint.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    pins['../'+NOTE.name] = sha256(NOTE.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


class MidpointBoxes:
    """O(r) exact preparation; O(1) integer operations per rectangle.

    Here r=min(weight,2*pairs-weight). Complementing a support subtracts
    weight-pairs from both full-pair counts when weight>pairs.
    """

    def __init__(self, pairs, weight, *, counts=None, budget=None):
        budget = budget or midpoint.near.provider.old.Budget(60.)
        counts = counts or midpoint.MidpointCounts(pairs, weight, budget=budget)
        if type(counts) is not midpoint.MidpointCounts or (counts.pairs, counts.weight) != (pairs, weight):
            raise ValueError('Expected matching midpoint count object')
        self.pairs, self.weight, self.counts = pairs, weight, counts
        self.r, self.shift, self.last = counts.r, counts.shift, counts.low
        self.denominator = counts.denominator
        self.v0, self.v1, self.e, self.f, self.triangle, self.diagonal = ([] for _ in range(6))
        for prefix in (self.v0, self.v1, self.e, self.f, self.triangle, self.diagonal):
            prefix.append(0)
        low, high = counts.low, counts.high
        constant, linear = low*high*(pairs-self.r), low*high
        for a, (left, right) in enumerate(zip(counts.left, counts.right)):
            if not a % 64:
                budget.check()
            self.v0.append(self.v0[-1]+right)
            self.v1.append(self.v1[-1]+a*right)
            scale = left*(self.r-2*a+1)
            e, f = scale*(constant+linear*a), scale*(linear-pairs*a)
            self.e.append(self.e[-1]+e)
            self.f.append(self.f[-1]+f)
            self.triangle.append(self.triangle[-1]+e*self.v0[-1]+f*self.v1[-1])
            self.diagonal.append(self.diagonal[-1]+right*(e+a*f))
        # Make completed tables immutable before returning them to a caller.
        for name in ('v0', 'v1', 'e', 'f', 'triangle', 'diagonal'):
            setattr(self, name, tuple(getattr(self, name)))
        budget.check()

    def _triangle_prefix(self, x, y):
        """Unscaled numerator over 0<=a<=x, 0<=b<=min(a,y)."""
        x, y = min(x, self.last), min(y, self.last)
        if min(x, y) < 0:
            return 0
        if x <= y:
            return self.triangle[x+1]
        return (self.triangle[y+1]
                + self.v0[y+1]*(self.e[x+1]-self.e[y+1])
                + self.v1[y+1]*(self.f[x+1]-self.f[y+1]))

    def _triangle_rectangle(self, al, ah, bl, bh):
        prefix = self._triangle_prefix
        return prefix(ah, bh)-prefix(al-1, bh)-prefix(ah, bl-1)+prefix(al-1, bl-1)

    def rectangle(self, al, ah, bl, bh):
        """Sum midpoint_upper(a,b) on an inclusive original-coordinate box.

        Out-of-domain coordinates contribute zero. An empty interval has
        sum zero. The diagonal is counted once, including its exceptional
        alternating support when r=pairs and a=b=0 after complementation.
        """
        if any(type(value) is not int for value in (al, ah, bl, bh)):
            raise ValueError('Integer rectangle endpoints required')
        al, ah = max(0, al-self.shift), min(self.last, ah-self.shift)
        bl, bh = max(0, bl-self.shift), min(self.last, bh-self.shift)
        if al > ah or bl > bh:
            return 0
        if self.r <= 1:
            return 1 if self.r == 0 else 2*self.pairs
        # The numerator is symmetric in a,b: left[a]*right[b] equals
        # left[b]*right[a], also for odd r (their array ratio is constant).
        numerator = (self._triangle_rectangle(al, ah, bl, bh)
                     + self._triangle_rectangle(bl, bh, al, ah))
        dl, dh = max(al, bl), min(ah, bh)
        if dl <= dh:
            numerator -= self.diagonal[dh+1]-self.diagonal[dl]
        value, remainder = divmod(self.pairs*numerator, self.denominator)
        if remainder or value < 0:
            raise ArithmeticError('Invalid midpoint rectangle integer sum')
        if self.r == self.pairs and al == bl == 0:
            value += 2
        return value

    def diagnostics(self):
        return dict(method='exact rectangle sum of raw midpoint upper counts',
                    pairs=self.pairs, active_groups=self.weight, canonical_weight=self.r,
                    profile_pair_table_built=False, pointwise_minimum_summed=False,
                    integer_prefix_arrays=6, integer_prefix_entries=6*(self.last+2),
                    rectangle_query_integer_operations='constant; bigint operand sizes depend on pairs',
                    source_sha256=source_pins())
