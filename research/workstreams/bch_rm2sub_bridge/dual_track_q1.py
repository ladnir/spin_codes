"""Matched BCH-256 diagnostics: transfer improvement versus spectrum uncertainty.

Uses the certified t64/s20 map, not the estimator's different nested map.
Every reported margin is Q1-only and binary64 until separately certified.
"""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
import sys

import numpy as np
from scipy.special import logsumexp

import bridge as base
import larger_state_maps as maps

LANDSCAPE = base.ROOT / 'workstreams/finite_asymptotic_theory/landscape_db'
sys.path.insert(0, str(LANDSCAPE))
import activation_q1 as old
import activation_q1_refresh as refresh
import bch_growth_model as onset


def even_binomial():
    weights = list(range(38, 220, 2))
    mass = sum(math.comb(256, w) for w in weights)
    counts = {w: F(((1 << 128) - 2) * math.comb(256, w), mass) for w in weights}
    return {0: F(1), **counts, 256: F(1)}


def aggregate(counts, coefficients):
    weights = [w for w, n in counts.items() if w and n]
    terms = np.array([math.log(counts[w]) + coefficients[w] for w in weights])
    return dict(margin_bits=-float(logsumexp(terms)) / math.log(2),
                dominant_weight=weights[int(np.argmax(terms))])


def values(method, spectrum, rows, log_a):
    lam = np.exp(log_a) / rows
    if method == 'refresh':
        moment = refresh.coefficient_logs(*refresh.epoch_logs(64, 20, spectrum, lam), rows // 64, 256)
    else:
        moment = old.coefficient_logs(*old.region_logs(*old.epoch_logs(64, 20, spectrum, lam), rows // 64), 256)
    return math.log(rows) + np.minimum(0., moment + (256 * rows // 10) * lam[:, None])


def run(exponents, output):
    assert not output.exists(), 'Use a fresh diagnostic output path'
    _, _, spectrum, _ = maps.load('t64_s20')
    counts = even_binomial()
    floor = {38: 1095168, 40: 18045952, 42: 42289664, 44: 267281408,
             46: 1247821056, 48: 5525363968, 50: 1574477312}
    floor.update({256 - w: n for w, n in list(floor.items())})
    floor[256] = 1
    model = onset.spectrum_model(256, counts)
    for shell in model['shells']:
        shell['count'] = float(shell['count'])
        if not math.isfinite(shell['log_upper']):
            shell['log_upper'] = None
    results = []
    for exponent in exponents:
        rows = (1 << exponent) // 128
        if rows < 64 or rows % 64:
            raise ValueError('Whole 64-bit epochs required')
        grid = np.arange(-40, 81) / 10
        coarse = values('refresh', spectrum, rows, grid)
        chosen = np.argmin(coarse, axis=0)
        # Refine low shells used by the spectrum estimate and the rigorous cap
        # envelope; reuse the same witness set for the old/new comparison.
        fine = np.unique(np.concatenate([grid[chosen[w]] + np.arange(-30, 31) / 500
                                         for w in (38, 40, 42, 44, 46, 48, 50)]))
        refined = values('refresh', spectrum, rows, fine)
        all_grid = np.r_[grid, fine]
        all_values = np.concatenate((coarse, refined))
        index = np.argmin(all_values, axis=0)
        after = all_values[index, np.arange(257)]
        before = values('old', spectrum, rows, all_grid).min(axis=0)
        if np.any(after[38:] > before[38:] + 1e-7):
            raise ArithmeticError('Refresh is weaker on a matched shell')
        row = dict(message_exponent=exponent, message_bits=1 << exponent,
                   outer_rows=rows, cutoff=256 * rows // 10,
                   modeled=aggregate(counts, after), old_modeled=aggregate(counts, before),
                   discovered_words_only=aggregate(floor, after),
                   witness_log_a={str(w): float(all_grid[index[w]]) for w in base.WEIGHTS},
                   modeled_intercept_bits=aggregate(counts, after)['margin_bits'] + math.log2(rows),
                   onset_model_margin_bits=model['onset_intercept_bits'] - math.log2(rows),
                   onset_chernoff_margin_bits=model['chernoff_intercept_bits'] - math.log2(rows))
        for label, coefficients in [('old_caps', before), ('refresh_caps', after)]:
            rational = {w: F.from_float(math.exp(max(-700., float(coefficients[w])))) for w in base.WEIGHTS}
            upper, _, _ = base.bch_bound(rational)
            row[label + '_margin_bits'] = math.log2(upper.denominator) - math.log2(upper.numerator)
        row['stress'] = []
        for bits in (1, 4, 8, 12, 16):
            stressed = {w: n * ((1 << bits) if (38 <= w <= 50 or 206 <= w <= 218) else 1)
                        for w, n in counts.items()}
            row['stress'].append(dict(low_shell_inflation_bits=bits, **aggregate(stressed, after)))
        results.append(row)
        print(f"k=2^{exponent}: caps {row['old_caps_margin_bits']:.6f} -> {row['refresh_caps_margin_bits']:.6f}; "
              f"modeled {row['modeled']['margin_bits']:.6f}; onset {row['onset_model_margin_bits']:.6f}", flush=True)
    paths = [Path(__file__), Path(old.__file__), Path(refresh.__file__), Path(onset.__file__),
             base.HERE / 'generated/larger_state_inputs_v1/t64_s20_selection.json',
             base.BCH / 'SPIN_HEURISTIC_EVIDENCE.md']
    base.write_new(output, dict(status='MATCHED_BCH256_Q1_DUAL_TRACK_DIAGNOSTIC',
        configuration='t64_s20', rows=results, onset_model=model,
        source_sha256={p.relative_to(base.ROOT).as_posix(): base.sha(p) for p in paths},
        limitations=['Nearest binary64, Q1 only; no new certificate in this file.',
                     'BCH cap inputs are rigorous but their transfer evaluation here is diagnostic.',
                     'Conditioned-even binomial spectrum and inflated-tail scenarios are not known code spectra.',
                     'Discovered-words-only is a partial envelope sum, not a lower bound on actual setup failure.',
                     'Onset model omits cancellations and finite-epoch fluctuations; not an encoder failure bound.',
                     'Higher occupancies have not been recomputed at larger message sizes.']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponents', nargs='+', type=int, default=[16, 20, 24, 28, 32])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.exponents, args.output)
