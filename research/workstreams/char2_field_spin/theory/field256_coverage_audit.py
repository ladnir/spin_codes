"""K20 dense-class certificate and integer cutoff audit for GRS[256,176].

This does not execute the separate sparse or intermediate replays. Rational
tilts are proposed numerically; only exact rational inequalities accept them.
No native code, timing, or extension-field word enumeration is used.
"""

import argparse
from fractions import Fraction
import json

from cauchy_firstmoment_screen import geometry
from regional_global_route import (GRID, choose_tilt,
                                   certify_density_interval, interval_base)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def audit(show_tilts=False):
    g = geometry(20, 256, 176)
    rows, n = g['rows'], g['n']
    require((rows, g['k'], n, g['threshold']) == (6016, 1058816, 1540096, 154009),
            'unexpected finite geometry')
    middle_first = (3*n+99)//100
    dense_first = (6*n+24)//25
    next_active_minimum = 81*601
    require(middle_first == 46203 and dense_first == 369624,'unexpected weight cutoffs')
    require(next_active_minimum == 48681 and next_active_minimum >= middle_first,
            'sparse and middle classes leave a gap')

    tilts = []
    for index in range(48, GRID):
        tilt = choose_tilt(32, 22, index)
        require(certify_density_interval(32, 22, index, tilt),
                f'failed small-outer interval {index}')
        lower, upper = Fraction(index, GRID), Fraction(index+1, GRID)
        small = interval_base(32, 22, lower, upper, tilt)
        large = interval_base(256, 176, lower, upper, tilt)
        # The large bracket is X**8+Y**8; the small bracket is X+Y.
        require(large <= small**8,f'failed eighth-power comparison {index}')
        require(certify_density_interval(256, 176, index, tilt),
                f'failed large-outer interval {index}')
        tilts.append((index, tilt.numerator, tilt.denominator))

    # The q=256 radius-.1 ball is bounded by 2**(1.27*N).
    require(Fraction(10, 9)**100 * 2295**10 < 2**127,'failed output-ball bound')
    # rows+1 <= 2**rows.bit_length(); N/100 >= floor(N/100).
    require(rows+1 <= 2**rows.bit_length(),'failed row-conditioning bound')
    margin = n//100-256*rows.bit_length()
    require(margin == 12072,'unexpected dense margin')
    result = dict(outer=[256, 176], q=256, width=16, state=8, k_log=20,
                  rows=rows, padded_source=g['k'], output_length=n,
                  forbidden_output_threshold=g['threshold'],
                  sparse_active_rows=[1, 600], sparse_replay_checked=False,
                  intermediate_weights=[middle_first, dense_first-1],
                  intermediate_replay_checked=False,
                  minimum_weight_if_active_rows_ge_601=next_active_minimum,
                  dense_weights=[dense_first, n],
                  dense_intervals_checked=len(tilts),
                  membership_exponent='32/25', ball_exponent='127/100',
                  dense_margin_bits_lower=margin, dense_rigorous=True,
                  cutoff_cover_complete=True, whole_code_certified=False)
    if show_tilts:
        result['rational_tilts'] = tilts
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--show-tilts', action='store_true')
    args = parser.parse_args()
    print(json.dumps(audit(args.show_tilts)), flush=True)
