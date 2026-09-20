"""Retain bounded syndrome density when an IMT state is activated.

If X has weight j and q=0, each target state has weighted mass at most
cap_j * z^j / C(t,j). Charge that amount to each nonzero target instead of
forgetting the distribution in D. The two representations are alternatives;
we select one complete zero row for each fixed j, never take entrywise minima.
"""
from functools import lru_cache
from fractions import Fraction as F
import math

import ladder

model = ladder.model


def zero_row(spectrum, kernel_count, total, cap, moment):
    levels = sorted(spectrum)
    return tuple([model.up(model.number(F(kernel_count, total)) * moment),
                  model.number(F(0))] +
                 [model.up(model.number(F(spectrum[w] * cap, total)) * moment)
                  for w in levels])


class Engine(ladder.Engine):
    @lru_cache(maxsize=24)
    def epoch(self, log_lam):
        rows = super().epoch(log_lam)
        z = (-model.number(log_lam).exp()).exp()
        output = []
        for j, old in enumerate(rows):
            total, cap = math.comb(128, j), int(self.caps[j]['cap'])
            if self.m * cap <= 4 * total and total > self.kernel[j]:
                row = zero_row(self.spectrum, self.kernel[j], total, cap, z**j)
                output.append(row + old[self.n:])
            else:
                output.append(old)
        return output
