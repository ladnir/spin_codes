"""Reuse fixed swap alias atlases across a complete interval of weights.

Fix z, alpha, a finite bank of positive coefficient markers, and one
positive state potential v. Each atlas is independent of the target
histogram. At weight r and full-pair count f, normalize its selected
residue by the individual multinomial and marker monomial. Only then
take entrywise minima across the bank.

The run-convolution bound sums every matching profile at each weight.
This helper evaluates every integer weight in the requested interval;
it does not interpolate from samples. All arithmetic here is FLOAT.
The FFT cushion is not a directed error bound or a certificate.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import isfinite, log2
from pathlib import Path

import numpy as np

from hd_grouped_packet_screen import LIMIT
from hd_shifted_bulk_search import LOG_DENSITY, LOG_FACTORIAL
from hd_shifted_projected_screen import alias_table
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
    checked_outer_measure, log_density_float, bound_description)
from hd_shifted_run_potential import run_log_screen
from hd_shifted_swap_potential import full_counts
from hd_shifted_swap_screen import CONSTRUCTION


class FixedAtlasScreen:
    """Build each distinct atlas once; never refit markers or the potential."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, grid=128, _tables=None,
                 outer_measure=DEFAULT_OUTER_MEASURE):
        self.outer_measure = checked_outer_measure(outer_measure)
        self.log_density = log_density_float(self.outer_measure)
        self.z, self.alpha, self.log_ratio = map(float, (z, alpha, log_ratio))
        self.odds = tuple(dict.fromkeys(tuple(map(float, row)) for row in marker_log_odds))
        if (not 0 < self.z < 1 or not 0 < self.alpha <= 1 or not isfinite(self.log_ratio)
                or not isinstance(grid, int) or grid < 64 or grid & (grid-1)
                or not self.odds or len(self.odds) > 100
                or any(len(row) != 2 or not all(map(isfinite, row)) for row in self.odds)):
            raise ValueError('A fixed valid witness and at most100 finite marker pairs are required')
        self.grid = grid
        self.tables = tuple(alias_table(self.z, self.alpha, row, grid) if _tables is None
                            else _tables[row] for row in self.odds)
        # Negative infinity represents an exact zero in an injected log table.
        if any(table.shape != (4, grid, grid) or np.any(np.isnan(table))
               or np.any(np.isposinf(table)) for table in self.tables):
            raise ArithmeticError('Invalid fixed alias table')
        self.constant = -self.alpha*(log2(255)+LIMIT*log2(self.z))
        self.source_witness = None
        self.central_marker_fits = 0

    @classmethod
    def from_witness(cls, witness, *, grid=None, outer_measure=None):
        """Reconstruct a replayable fixed bank from exact literals, not convenience floats."""
        if (witness.get('construction') != CONSTRUCTION or witness.get('pair_swaps') is not True
                or witness.get('orientation_aggregation') != 'independent-fair-swaps'):
            raise ValueError('The witness must explicitly select independent fair swaps')
        saved_measure = checked_outer_measure(witness.get('outer_measure', DEFAULT_OUTER_MEASURE))
        measure = saved_measure if outer_measure is None else checked_outer_measure(outer_measure)
        source = witness.get('source_exact_witness')
        if (not isinstance(source, dict)
                or any(not isinstance(source.get(name), str) for name in ('z', 'alpha', 'potential'))):
            raise ValueError('Exact saved witness literals are required; raw FLOAT witnesses are not replayable')
        try:
            z, alpha, potential = (F(source[name]) for name in ('z', 'alpha', 'potential'))
            markers = tuple(tuple(F(value) for value in row.split()) for row in source['markers'])
            if (not 0 < z < 1 or not 0 < alpha <= 1 or potential <= 0
                    or not markers or len(markers) > 100
                    or any(len(row) != 2 or min(row) <= 0 for row in markers)):
                raise ValueError('Admissible exact fixed-witness literals are required')
            odds = tuple(tuple(log2(value) for value in row) for row in markers)
            ratio = log2(potential)
        except (KeyError, TypeError, ZeroDivisionError, OverflowError, AttributeError) as error:
            raise ValueError('Malformed exact fixed witness') from error
        result = cls(z, alpha, odds, ratio, grid=witness.get('grid', 128) if grid is None else grid,
                     outer_measure=measure)
        result.source_witness = deepcopy(source)
        return result

    @classmethod
    def from_checkpoint(cls, weight, *, grid=128, outer_measure=DEFAULT_OUTER_MEASURE):
        """Load exact saved literals, with no central or neighboring fit."""
        from hd_shifted_swap_potential import SAVED_RUN_WITNESSES
        if weight not in SAVED_RUN_WITNESSES:
            raise ValueError('No saved fixed atlas witness at this weight')
        source = SAVED_RUN_WITNESSES[weight]
        odds = tuple(tuple(log2(F(value)) for value in row.split()) for row in source['markers'])
        result = cls(F(source['z']), F(source['alpha']), odds,
                     log2(F(source['potential'])), grid=grid, outer_measure=outer_measure)
        result.source_witness = dict(source, active_groups=weight)
        return result

    @classmethod
    def fit_central_once(cls, weight, z, alpha, centers, log_ratio, *, grid=128,
                         outer_measure=DEFAULT_OUTER_MEASURE):
        """Fit a finite bank at one weight, then freeze it for all evaluations.

        Each center uses the established atlas_bank_screen proposer once.
        Retain each generated table before the shared LRU cache can evict
        it; neighboring weights never rebuild these tables or fit markers.
        """
        from hd_shifted_swap_potential import atlas_bank_screen
        outer_measure = checked_outer_measure(outer_measure)
        centers = tuple(sorted(set(centers)))
        if (not isinstance(weight, int) or not 1 <= weight <= 8192
                or not 0 < z < 1 or not 0 < alpha <= 1 or not isfinite(log_ratio)
                or not isinstance(grid, int) or grid < 64 or grid & (grid-1)
                or not centers or len(centers) > 100
                or any(f not in full_counts(4096, weight) for f in centers)):
            raise ValueError('At most100 feasible explicit central marker targets are required')
        z, alpha, log_ratio = map(float, (z, alpha, log_ratio))
        odds, tables, evaluations = [], {}, 0
        for full in centers:
            _, fitted = atlas_bank_screen(weight, z, alpha, (full,), grid=grid)
            row = tuple(fitted[0]['odds'])
            odds.append(row)
            # This query hits the table just generated by atlas_bank_screen.
            tables[row] = alias_table(z, alpha, row, grid)
            evaluations += fitted[0]['objective_evaluations']
        result = cls(z, alpha, odds, log_ratio, grid=grid, _tables=tables, outer_measure=outer_measure)
        result.central_marker_fits = len(centers)
        result.source_witness = dict(
            active_groups=weight, z=str(F(z)), alpha=str(F(alpha)),
            potential=str(F(2**log_ratio)),
            markers=[' '.join(str(F(2**value)) for value in row) for row in result.odds],
            central_targets=list(centers), central_objective_evaluations=evaluations,
            marker_exactification='exact-rational-value-of-binary64-exp2')
        return result

    def regional(self, weight):
        """Return normalized log2 regional bounds for every feasible full count."""
        counts = tuple(full_counts(4096, weight))
        full = np.asarray(counts, dtype=np.int64)
        single, empty = weight-2*full, 4096-weight+full
        log_mult = (LOG_FACTORIAL[4096]-LOG_FACTORIAL[empty]
                    -LOG_FACTORIAL[single]-LOG_FACTORIAL[full])
        best = np.full((4, len(full)), np.inf)
        for row, table in zip(self.odds, self.tables):
            normalized = (table[:, single % self.grid, full % self.grid]
                          -log_mult-single*row[0]-full*row[1])
            np.minimum(best, normalized, out=best)
        if np.any(np.isnan(best)) or np.any(np.isposinf(best)):
            raise ArithmeticError('Nonfinite normalized regional screen')
        return {int(f): best[:, index].reshape(2, 2) for index, f in enumerate(full)}

    def evaluate(self, weight):
        """Sum all indexed scalar supports at one nonzero active-group weight."""
        if not isinstance(weight, int) or not 1 <= weight <= 8192:
            raise ValueError('A nonzero integer scalar support weight is required')
        regions = self.regional(weight)
        support_sum = run_log_screen(4096, weight, regions, self.log_ratio)
        prefix = self.alpha*self.log_density*weight+self.constant
        return dict(active_groups=weight, regional_histograms=len(regions),
                    log2_support_sum_screen=support_sum,
                    log2_outer_prefix=prefix, log2_upper_screen=prefix+support_sum,
                    covers_all_matching_profiles=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False)

    def witness(self):
        result = dict(mark=self.z, alpha=self.alpha, marker_log_odds=self.odds,
                    log2_potential_ratio=self.log_ratio, grid=self.grid,
                    construction=CONSTRUCTION, pair_swaps=True,
                    orientation_aggregation='independent-fair-swaps',
                    regional_minimum='after-histogram-normalization',
                    source_exact_witness=self.source_witness)
        if self.outer_measure != DEFAULT_OUTER_MEASURE:
            result.update(outer_measure=self.outer_measure,
                          outer_measure_bound=bound_description(self.outer_measure))
        return result

    def band(self, first, last, *, threshold=-40., max_weights=33, progress=None):
        """Evaluate the entire inclusive interval, subject to an explicit cap."""
        if (not isinstance(first, int) or not isinstance(last, int)
                or not 1 <= first <= last <= 8192
                or not isinstance(max_weights, int) or max_weights < 1
                or last-first+1 > max_weights or not isfinite(threshold)):
            raise ValueError('A nonempty integer interval within the explicit weight cap is required')
        results = []
        for weight in range(first, last+1):
            row = self.evaluate(weight)
            results.append(row)
            if progress:
                progress(dict(stage='fixed-atlas-weight', **row))
        failed = [row['active_groups'] for row in results if row['log2_upper_screen'] > threshold]
        union = float(np.logaddexp2.reduce([row['log2_upper_screen'] for row in results]))
        result = dict(requested_weight_range=[first, last], weights_evaluated=len(results),
                    max_weights=max_weights, every_weight_evaluated=True,
                    all_weights_pass_threshold=not failed, failed_weights=failed,
                    threshold=threshold, union_log2_upper_screen=union,
                    union_passes_threshold=union <= threshold, results=results,
                    fixed_witness=self.witness(), fft_tables_built=len(self.tables),
                    central_marker_fits=self.central_marker_fits,
                    witness_searches=self.central_marker_fits, neighbor_witness_searches=0,
                    directed=False, whole_code_certificate=False,
                    construction=CONSTRUCTION, covers_all_matching_profiles=True,
                    covers_all_even_odd_splits=True,
                    method='fixed-atlas-common-potential-run-convolution')
        if self.outer_measure != DEFAULT_OUTER_MEASURE:
            result.update(outer_measure=self.outer_measure,
                          outer_measure_bound=bound_description(self.outer_measure))
        return result


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--witness', type=Path,
                           help='Frozen JSON with mark, alpha, marker_log_odds and log2_potential_ratio')
    selection.add_argument('--checkpoint', type=int, help='Use exact SAVED_RUN_WITNESSES literals')
    selection.add_argument('--fit-frontier-1024', action='store_true',
                           help='Fit the stronger fixed1024 frontier bank once, never neighboring weights')
    parser.add_argument('--first', type=int, default=1008)
    parser.add_argument('--last', type=int, default=1040)
    parser.add_argument('--max-weights', type=int, default=33)
    parser.add_argument('--threshold', type=float, default=-40.)
    parser.add_argument('--grid', type=int, help='Override the alias grid; default128 for checkpoints')
    parser.add_argument('--outer-measure', choices=(DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE),
                        help='Select a proved outer bound; defaults to the saved mode or inherited-diagonal')
    parser.add_argument('--output', type=Path, help='Exclusive-create a new FLOAT result artifact')
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    if args.fit_frontier_1024:
        centers = tuple(sorted(set(range(0, 513, 32)) | {512-k for k in (1, 2, 4, 8, 16)}))
        screen = FixedAtlasScreen.fit_central_once(1024, .4828132651072611, .10118724647150915,
            centers, -8.532552444663622, grid=128 if args.grid is None else args.grid,
            outer_measure=args.outer_measure or DEFAULT_OUTER_MEASURE)
    elif args.checkpoint is not None:
        screen = FixedAtlasScreen.from_checkpoint(args.checkpoint,
            grid=128 if args.grid is None else args.grid,
            outer_measure=args.outer_measure or DEFAULT_OUTER_MEASURE)
    else:
        with args.witness.open(encoding='utf-8') as source:
            witness = json.load(source)
        screen = FixedAtlasScreen.from_witness(witness, grid=args.grid, outer_measure=args.outer_measure)
    result = screen.band(args.first, args.last, threshold=args.threshold, max_weights=args.max_weights,
                         progress=lambda row: print(json.dumps(row), flush=True))
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            json.dump(result, target, separators=(',', ':'), allow_nan=False)
            target.write('\n')
        print(json.dumps(dict(path=str(args.output.resolve()),
                              **{key: value for key, value in result.items()
                                 if key not in ('results', 'fixed_witness')})), flush=True)
    else:
        print(json.dumps(result, allow_nan=False), flush=True)
