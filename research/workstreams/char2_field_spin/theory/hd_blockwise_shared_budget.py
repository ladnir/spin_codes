"""All-route quiet-subcode cap for tile-local invertible feedthrough.

This is a necessary design gate, not a lower distance certificate. The
averaging is over a choice of source rows after fixing the route; no random
route assumption enters the theorem.
"""

from fractions import Fraction
from math import comb


def covered_mean(rows, active, tile, output):
    if tile < 1 or rows < 1 or rows % tile or not 0 <= active <= rows or output % rows:
        raise ValueError('whole regional tiles and a valid row subset required')
    miss = Fraction(comb(rows-active, tile), comb(rows, tile))
    return output*(1-miss)


def tiled_quiet_cap(message_symbols=1 << 20, outer_length=256,
                    outer_dimension=128, width=64, state=4, tile=8):
    if (message_symbols < 1 or not 1 <= outer_dimension <= outer_length
            or message_symbols % outer_dimension or width < 1 or state < 0
            or tile < 1 or width % tile):
        raise ValueError('invalid dimensions')
    rows = message_symbols//outer_dimension
    output = rows*outer_length
    if rows % tile or output % width:
        raise ValueError('regional tiles and packets must be complete')
    equations = state*(output//width-1)
    active = equations//outer_dimension+1
    singleton = output-message_symbols+1
    result = dict(rows=rows, output_symbols=output, outer=[outer_length, outer_dimension],
                  tile=tile, width=width, state=state, selected_rows=active,
                  feedback_equations=equations, deterministic_every_route=True,
                  whole_certificate=False, quiet_witness=active <= rows)
    if active > rows:
        return dict(result, quiet_dimension=None, covered_output_mean=None,
                    quiet_weight_cap=None, distance_upper=str(Fraction(singleton, output)),
                    rejects_40_percent=Fraction(singleton, output) < Fraction(2, 5))
    dimension = active*outer_dimension-equations
    mean = covered_mean(rows, active, tile, output)
    cap = mean.numerator//mean.denominator-dimension+1
    if not 1 <= cap <= output:
        raise ArithmeticError('invalid quiet-subcode weight cap')
    return dict(result, quiet_dimension=dimension, covered_output_mean=str(mean),
                quiet_weight_cap=cap, distance_upper=str(Fraction(min(cap, singleton), output)),
                rejects_40_percent=Fraction(min(cap, singleton), output) < Fraction(2, 5))


def receipt():
    return [tiled_quiet_cap(width=t, state=s, tile=g)
            for t, s, g in ((64, 4, 1), (64, 1, 8), (64, 2, 8),
                            (64, 4, 4), (64, 4, 8), (128, 4, 8),
                            (256, 4, 8))]


if __name__ == '__main__':
    import json
    print(json.dumps(receipt(), indent=2))
