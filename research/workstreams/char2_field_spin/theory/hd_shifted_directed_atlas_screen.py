"""FLOAT proposals whose alias tables come from directed regional transforms.

The input is a fixed swap witness with exact rational literals. Build each
distinct marker atlas once, then reuse its logarithms at every requested
weight. The conversion restores each stored scaling factor in log space;
it neither inserts a numerical floor nor converts a huge factor to float.

Only the source transforms use directed arithmetic. Histogram normalization,
potential contractions, support sums, and the final union remain FLOAT.
Every result requires a fresh directed replay before certificate use.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import isfinite

import numpy as np
from flint import arb, ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_coverage import RegionalAliases
from hd_shifted_good_pair_certificate import component_power
from hd_shifted_outer_density import DEFAULT_OUTER_MEASURE, checked_outer_measure
from hd_shifted_potential_band_screen import FixedAtlasScreen
from hd_shifted_run_count_screen import run_log_screen_float
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION, swap_unmarked
from hd_shifted_unoriented_certificate import unsigned_components


TABLE_SOURCE = 'directed-regional-alias-upper-endpoints'
COUNT_METHOD = 'FLOAT-log-factorial-cyclic-run-counts'


def checked_window(first, last, *, max_weights=33, threshold=-40.):
    """Check the explicit screen budget before any directed atlas is built."""
    if (type(first) is not int or type(last) is not int
            or not 1 <= first <= last <= 8192
            or type(max_weights) is not int or not 1 <= max_weights <= 8192
            or last-first+1 > max_weights or not isfinite(threshold)):
        raise ValueError('A nonempty integer interval within the explicit weight cap is required')


def atlas_log_table(atlas, *, precision=256):
    """Convert entry 2*i+j to log2(upper)+log2(factor), without rescaling.

    RegionalAliases already stores rho**128 * vector[i]/vector[j]. Its
    factor is used exactly once. A zero upper endpoint remains -infinity.
    The returned binary64 logarithms are proposal values, not certificates.
    """
    upper = np.asarray(atlas.upper)
    if (type(precision) is not int or precision < 128
            or upper.ndim != 3 or upper.shape[0] != 4
            or upper.shape[1] != upper.shape[2] or upper.shape[1] < 1
            or not np.all(np.isfinite(upper)) or np.any(upper < 0)
            or len(atlas.factors) != 4):
        raise ValueError('Finite nonnegative two-dimensional atlas uppers are required')
    with ctx.workprec(precision):
        if any(not factor.is_finite() or not factor > 0 for factor in atlas.factors):
            raise ArithmeticError('Every atlas scaling factor must be certainly positive and finite')
        log_two = arb(2).log()
        offsets = np.asarray([float((factor.log()/log_two).upper()) for factor in atlas.factors])
        if not np.all(np.isfinite(offsets)):
            raise ArithmeticError('Atlas factor logarithm is outside the FLOAT proposal range')
        with np.errstate(divide='ignore', invalid='raise', over='raise'):
            table = np.log2(upper) + offsets[:, None, None]
    if np.any(np.isnan(table)) or np.any(np.isposinf(table)):
        raise ArithmeticError('Invalid alias logarithm')
    table.setflags(write=False)
    return table


def _checked_source(witness, grid, precision, bits, outer_measure):
    if (not isinstance(witness, dict) or witness.get('construction') != SWAP_CONSTRUCTION
            or witness.get('pair_swaps') is not True
            or witness.get('orientation_aggregation') != 'independent-fair-swaps'):
        raise ValueError('An explicit independent-fair-swap witness is required')
    saved_measure = checked_outer_measure(witness.get('outer_measure', DEFAULT_OUTER_MEASURE))
    measure = saved_measure if outer_measure is None else checked_outer_measure(outer_measure)
    grid = witness.get('grid', 128) if grid is None else grid
    if (type(grid) is not int or grid < 64 or grid & (grid-1)
            or type(precision) is not int or precision < 128
            or type(bits) is not int or bits < 16):
        raise ValueError('An admissible grid, source precision, and packet precision are required')
    source = witness.get('source_exact_witness')
    if (not isinstance(source, dict)
            or any(not isinstance(source.get(name), str) for name in ('z', 'alpha', 'potential'))
            or not isinstance(source.get('markers'), (list, tuple))
            or not source['markers'] or len(source['markers']) > 100
            or any(not isinstance(row, str) for row in source['markers'])):
        raise ValueError('Exact saved witness literals are required')
    try:
        z, alpha, potential = (F(source[name]) for name in ('z', 'alpha', 'potential'))
        markers = tuple(dict.fromkeys(tuple(F(value) for value in row.split())
                                     for row in source['markers']))
        if (not 0 < z < 1 or not 0 < alpha <= 1 or potential <= 0
                or any(len(row) != 2 or min(row) <= 0 for row in markers)):
            raise ValueError('Inadmissible exact witness')
        with ctx.workprec(precision):
            log_two = arb(2).log()
            odds = tuple(tuple(float(rational(value).log()/log_two) for value in row)
                         for row in markers)
            log_ratio = float(rational(potential).log()/log_two)
        if (not 0 < float(z) < 1 or not 0 < float(alpha) <= 1
                or not isfinite(log_ratio) or any(not all(map(isfinite, row)) for row in odds)):
            raise ValueError('Exact witness is outside the FLOAT proposal range')
        # The base class keys tables by binary64 log markers. Do not silently
        # merge distinct exact rows whose logarithms round to the same key.
        if len(set(odds)) != len(markers):
            raise ValueError('Distinct exact marker pairs collide in FLOAT log coordinates')
    except (TypeError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('Malformed exact witness') from error
    return deepcopy(source), z, alpha, markers, odds, log_ratio, grid, measure


class DirectedAtlasScreen(FixedAtlasScreen):
    """A fixed FLOAT screen with directed source tables, never a receipt."""

    @classmethod
    def from_witness(cls, witness, *, grid=None, precision=256, bits=60,
                     outer_measure=None, progress=None):
        """Build from a fixed_witness object; do not fit or alter its literals.

        Precision and bits describe the source tables only. The independent
        directed replayer chooses its own precision and packet cap settings.
        """
        source, z, alpha, markers, odds, ratio, grid, measure = _checked_source(
            witness, grid, precision, bits, outer_measure)
        tables = {}
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='directed-screen-packet-coefficients', precision=precision,
                              packet_mantissa_bits=bits))
            single, paired = swap_unmarked(z, alpha, precision, bits)
            for index, (markers_row, odds_row) in enumerate(zip(markers, odds)):
                pieces = unsigned_components(single, paired, markers_row)
                powered, vector, rho = component_power(*pieces, grid, dimensions=2)
                atlas = RegionalAliases.from_power(powered, vector, rho)
                tables[odds_row] = atlas_log_table(atlas, precision=precision)
                del powered, atlas, pieces
                if progress:
                    progress(dict(stage='directed-screen-atlas', precision=precision,
                                  completed=index+1, total=len(markers), grid=grid))
        result = cls(z, alpha, odds, ratio, grid=grid, _tables=tables, outer_measure=measure)
        result.source_witness = source
        result.table_source_metadata = dict(kind=TABLE_SOURCE, precision=precision,
            packet_mantissa_bits=bits, distinct_markers=len(markers),
            factor_conversion='log2-upper-plus-arb-log2-factor', zero_alias='negative-infinity')
        return result

    def witness(self):
        result = super().witness()
        result.update(source_exact_witness=deepcopy(self.source_witness),
                      outer_measure=self.outer_measure,
                      table_source=deepcopy(self.table_source_metadata),
                      count_method=COUNT_METHOD,
                      requires_fresh_directed_replay=True)
        return result

    def evaluate(self, weight):
        """Use the existing FLOAT log-factorial counts, without bigint tables."""
        if type(weight) is not int or not 1 <= weight <= 8192:
            raise ValueError('A nonzero integer scalar support weight is required')
        regions = self.regional(weight)
        support_sum = run_log_screen_float(4096, weight, regions, self.log_ratio)
        prefix = self.alpha*self.log_density*weight+self.constant
        return dict(active_groups=weight, regional_histograms=len(regions),
                    log2_support_sum_screen=support_sum, log2_outer_prefix=prefix,
                    log2_upper_screen=prefix+support_sum,
                    covers_all_matching_profiles=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False,
                    requires_fresh_directed_replay=True, table_source=TABLE_SOURCE,
                    count_method=COUNT_METHOD)

    def band(self, first, last, **kwargs):
        checked_window(first, last, max_weights=kwargs.get('max_weights', 33),
                       threshold=kwargs.get('threshold', -40.))
        return dict(super().band(first, last, **kwargs),
                    requires_fresh_directed_replay=True,
                    outer_measure=self.outer_measure,
                    count_method=COUNT_METHOD,
                    table_source=deepcopy(self.table_source_metadata))


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--witness', type=Path, required=True,
                        help='Exact fixed_witness JSON, or a proposal containing fixed_witness')
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--max-weights', type=int, default=33)
    parser.add_argument('--threshold', type=float, default=-40.)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=60)
    parser.add_argument('--grid', type=int)
    parser.add_argument('--outer-measure')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        checked_window(args.first, args.last, max_weights=args.max_weights, threshold=args.threshold)
    except ValueError as error:
        parser.error(str(error))
    if args.output and args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    with args.witness.open(encoding='utf-8') as source_file:
        witness = json.load(source_file)
    witness = witness.get('fixed_witness', witness)
    progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
    screen = DirectedAtlasScreen.from_witness(witness, precision=args.precision, bits=args.bits,
        grid=args.grid, outer_measure=args.outer_measure, progress=progress)
    result = screen.band(args.first, args.last, max_weights=args.max_weights,
                         threshold=args.threshold, progress=progress)
    # Fail closed on non-JSON numeric output. Genuine zero table entries are
    # retained internally; no floor is introduced to make a result serializable.
    encoded = json.dumps(result, allow_nan=False, separators=(',', ':'))
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            target.write(encoded+'\n')
        print(json.dumps(dict(path=str(args.output.resolve()),
                              union_log2_upper_screen=result['union_log2_upper_screen'],
                              failed_weights=result['failed_weights'], directed=False,
                              requires_fresh_directed_replay=True)), flush=True)
    else:
        print(encoded, flush=True)
