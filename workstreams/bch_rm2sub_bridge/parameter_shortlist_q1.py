"""Serial, read-only Q1 diagnostics for the BCH-256 parameter shortlist.

Not an outward certificate or an all-occupancy check. Requires retained maps
and BCH inputs; prints results without changing frozen files or receipts.
"""
import argparse
import math
from fractions import Fraction as F

import numpy as np

import bridge as base
import dual_track_q1 as diagnostic
import larger_state_maps as maps

CONFIGURATIONS = (*maps.NAMES, *base.CONFIGS)


def margin(t, s, spectrum, exponent):
    rows = (1 << exponent) // 128
    if exponent < 7 or rows < t or rows % t:
        raise ValueError('An integral number of complete epochs per region is required')
    grid = np.arange(-40, 81) / 10

    def evaluate(log_a):
        lam = np.exp(log_a) / rows
        moments = diagnostic.refresh.coefficient_logs(
            *diagnostic.refresh.epoch_logs(t, s, spectrum, lam), rows // t, 256)
        return math.log(rows) + np.minimum(
            0., moments + (256 * rows // 10) * lam[:, None])

    coarse = evaluate(grid)
    chosen = np.argmin(coarse, axis=0)
    fine = np.unique(np.concatenate([
        grid[chosen[w]] + np.arange(-30, 31) / 500
        for w in (38, 40, 42, 44, 46, 48, 50)]))
    coefficients = np.minimum(coarse.min(axis=0), evaluate(fine).min(axis=0))
    upper, _, _ = base.bch_bound({
        w: F.from_float(math.exp(max(-700., float(coefficients[w]))))
        for w in base.WEIGHTS})
    return math.log2(upper.denominator) - math.log2(upper.numerator)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--configurations', nargs='+', choices=CONFIGURATIONS,
                        default=['t64_s20', 't128_s19', 't64_s16'])
    parser.add_argument('--exponents', nargs='+', type=int, default=[16, 18, 20])
    args = parser.parse_args()
    for name in args.configurations:
        if name in maps.NAMES:
            t, s, spectrum, _ = maps.load(name)
        else:
            t, s, spectrum = base.load_map(name)
        for exponent in args.exponents:
            value = margin(t, s, spectrum, exponent)
            print(f'{name} K=2^{exponent} Q1_CAP_DIAGNOSTIC {value:.6f}', flush=True)
