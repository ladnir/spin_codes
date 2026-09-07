"""Bound activation-state density using the exact character weight spectrum.

J is the transpose of the recorded generator of A. Character inversion
bounds each syndrome probability by an absolute Krawtchouk sum. When
cheap enough, replace the arbitrary activation class by a dilated L class.
"""
import math

import numpy as np

import occupation_refresh_v1 as base


def krawtchouk_row(length, weight):
    values = [1, length-2*weight]
    for j in range(1, length):
        numerator = (length-2*weight)*values[j]-(length-j+1)*values[j-1]
        value, remainder = divmod(numerator, j+1)
        if remainder:
            raise ArithmeticError('Krawtchouk recurrence lost exactness')
        values.append(value)
    return values


class Epochs(base.Epochs):
    def __init__(self, t, s, a_counts, kernel, maximum_inflation=4):
        super().__init__(t, s, a_counts, kernel)
        field = 1 << s
        if sum(a_counts.values()) != field-1:
            raise ValueError('nonzero character spectrum has wrong total')
        rows = {w: krawtchouk_row(t, w) for w in a_counts}
        self.activation = []
        for j in range(t+1):
            choose = math.comb(t, j)
            signed = choose+sum(n*rows[w][j] for w, n in a_counts.items())
            if signed != field*kernel[j]:
                raise ValueError('kernel does not match the character spectrum')
            numerator = choose+sum(n*abs(rows[w][j]) for w, n in a_counts.items())
            active = choose-kernel[j]
            numerator = min(numerator, field*active)
            # L allows point probabilities <= 1/(2^s-2). Dilation must
            # dominate both total nonzero mass and each individual state.
            mass_numerator = max(field*active, (field-2)*numerator)
            use_density = j > 0 and active > 0 and mass_numerator <= maximum_inflation*field*active
            self.activation.append(dict(numerator=numerator, denominator=field*choose,
                                        mass_numerator=mass_numerator, use_density=use_density))

    def at(self, lam, maximum=None):
        result = super().at(lam, maximum).copy()
        for j, law in enumerate(self.activation[:len(result)]):
            if law['use_density']:
                result[j, 0, 1] = -np.inf
                result[j, 0, 3] = math.log(law['mass_numerator'])-math.log(law['denominator'])-lam*j
        return result
