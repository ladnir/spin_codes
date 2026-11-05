"""Selected-weight cell bounds from fully normalized regional cap tables.

The frozen RegionalCache normalizes every histogram and takes entrywise
marker minima before this module takes interval maxima. The selected weight
is fixed throughout; neither interpolation nor real-histogram envelopes are
used. A sparse table answers componentwise interval maxima in constant time.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_gamma_cell_certificate as numeric
import hd_field128_gamma_cell_counts as counts
import hd_field128_gamma_profile_certificate as profiles
import hd_field128_lift_alias_cover as restored
from hd_field128_gamma_cell_geometry import Cell


base, tail, provider = numeric.base, numeric.tail, numeric.provider


def source_pins():
    directory, pending, pins = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (directory/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in pins:
                    pending.append(target)
    return dict(sorted(pins.items()))


def _maximum(left, right):
    return (max(left[0], right[0]), max(left[1], right[1]),
            max(left[2], right[2]), max(left[3], right[3]))


class RegionalMaxTable:
    """Immutable O(n log n) table, with two-block O(1) interval queries."""

    def __init__(self, pairs, weight, regions, budget):
        full = tail.full_counts(pairs, weight)
        if not 1 <= weight < 2*pairs or set(regions) != set(full):
            raise ValueError('Complete regional range at one nonempty, non-full weight required')
        budget.check()
        rows = []
        for index, f in enumerate(full):
            if not index % 64:
                budget.check()
            row = tuple(base.upper(value) for value in regions[f])
            if len(row) != 4:
                raise ValueError('Each normalized regional matrix must have four entries')
            rows.append(row)
        levels, half, size = [tuple(rows)], 1, len(rows)
        while 2*half <= size:
            prior, following = levels[-1], []
            for start in range(size-2*half+1):
                if not start % 64:
                    budget.check()
                following.append(_maximum(prior[start], prior[start+half]))
            levels.append(tuple(following))
            half *= 2
        budget.check()
        self.pairs, self.weight, self.first = pairs, weight, full.start
        self.last, self.levels = full.stop-1, tuple(levels)

    def query(self, low, high):
        if (type(low) is not int or type(high) is not int
                or not self.first <= low <= high <= self.last):
            raise ValueError('Interval must lie inside the complete selected-weight table')
        degree = (high-low+1).bit_length()-1
        span = 1 << degree
        rows = self.levels[degree]
        return _maximum(rows[low-self.first], rows[high-span+1-self.first])


def full_intervals(cell, weight):
    """Exact parity projections at the locked weight; None if integer-empty."""
    if (not isinstance(cell, Cell) or type(weight) is not int
            or not 1 <= weight < 2*cell.pairs):
        raise ValueError('Expected a native cell at one nonempty, non-full weight')
    pairs = cell.pairs
    if cell.kind == 'interior':
        (rl, rh), (al, ah), (bl, bh) = cell.ranges
        if rl != weight or rh != weight:
            raise ValueError('Every interior cell must be locked to the table weight')
        low, high = max(0, weight-pairs), (weight-1)//2
        even, odd = (max(al, low), min(ah, high)), (max(bl, low), min(bh, high))
        return None if even[0] > even[1] or odd[0] > odd[1] else (even, odd)
    if cell.kind in ('even_face', 'odd_face'):
        (ul, uh), (kl, kh) = cell.ranges
        if weight % 2 or ul != weight//2 or uh != weight//2:
            raise ValueError('Every face cell must be locked to the even table weight')
        u = weight//2
        kl, kh = max(1, kl), min(kh, u, pairs-u)
        if kl > kh:
            return None
        pure, other = (u, u), (u-kh, u-kl)
        return (pure, other) if cell.kind == 'even_face' else (other, pure)
    raise ValueError('Unknown selected-weight cell kind')


def _cell_bound(table, cell, prefix, support_count, budget):
    """Synthetic-test arithmetic core; public receipts verify provider/count."""
    budget.check()
    if (not isinstance(table, RegionalMaxTable) or not isinstance(cell, Cell)
            or table.pairs != cell.pairs or type(support_count) is not int or support_count < 0):
        raise ValueError('Expected matching exact cell/table geometry and an integer count')
    intervals = full_intervals(cell, table.weight)
    if intervals is None:
        budget.check()
        return dict(zero=True, integer_empty=True, full_count_intervals=None,
                    log2_upper_rational=None, log2_upper_enclosure=None)
    if not support_count:
        raise ValueError('A nonempty native cell requires a positive support bound')
    left, right = (table.query(*interval) for interval in intervals)
    transfer = profiles.pair_row_sum128(left, right)
    if transfer.is_zero():
        result = dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None)
    else:
        value = tail.log2arb(arb(support_count))+prefix+tail.log2arb(transfer)
        result = dict(zero=False, log2_upper_rational=str(base.rational_ceiling(value)),
                      log2_upper_enclosure=str(value))
    budget.check()
    return dict(result, integer_empty=False, full_count_intervals=intervals,
                transfer=profiles.support_log(transfer),
                regional_log2_upper_rational=[[
                    None if value.is_zero() else str(base.rational_ceiling(tail.log2arb(value)))
                    for value in row] for row in (left, right)])


class Evaluator(numeric.Evaluator):
    """Frozen-cover-compatible evaluator locked to one selected weight.

    Construction computes the normalized table once. A failed construction
    raises; a failed cell returns an incomplete receipt with no cell bound.
    """

    def __init__(self, bank, weight, *, max_seconds=120.):
        if type(weight) is not int or not 1 <= weight < base.GROUPS:
            raise ValueError('Selected weight must lie in1..8191')
        budget = base.Budget(max_seconds)
        super().__init__(bank)
        self.pins = source_pins()
        with ctx.workprec(bank.precision):
            regions = restored.RegionalCache(bank).regional(weight, budget)
            table = RegionalMaxTable(base.PAIRS, weight, regions, budget)
            prefix = numeric.signed_upper(bank.constant+bank.slope*weight)
        budget.check()
        self.weight, self.table, self.prefix = weight, table, prefix
        self.preparation_seconds = monotonic()-budget.started
        self.check_pins()

    def check_pins(self):
        profiles.check_source_pins(self.pins)
        self.bank.check_pins()

    def receipt(self, cell, *, marker_pairs=None, max_seconds=60., failure_bits=64):
        if (not isinstance(cell, Cell) or cell.pairs != base.PAIRS
                or type(failure_bits) is not int or not 0 <= failure_bits <= 4096):
            raise ValueError('Expected a production native cell and valid probability budget')
        full_intervals(cell, self.weight)
        if marker_pairs is not None:
            raise ValueError('This table always uses every bank marker; marker-pair filters are unsupported')
        budget, bank = base.Budget(max_seconds), self.bank
        identity = dict(status='pending', cell=cell.descriptor(), construction=bank.construction,
                        state=bank.state, q=str(bank.q), target_distance=str(bank.distance),
                        forbidden_output_weight=bank.forbidden, precision=bank.precision,
                        witness=dict(z=str(bank.z), alpha=str(bank.alpha), grid=bank.grid, bits=bank.bits,
                                     markers=[[str(value) for value in marker] for marker in bank.markers]),
                        active_groups=self.weight, selected_weight_locked=True,
                        directed=False, complete_cell=False, whole_code_certificate=False,
                        log2_upper_rational=None, failure_bits=failure_bits,
                        cell_source_sha256=self.pins, regional_table_source_sha256=self.pins,
                        bank_source_sha256=bank.pins,
                        bound_method='normalized-regional-table-componentwise-interval-maxima',
                        normalization='frozen RegionalCache: exact histogram normalization before all-marker entrywise minimum',
                        regional_histograms=self.table.last-self.table.first+1,
                        regional_table_builds=1, regional_table_preparation_seconds=self.preparation_seconds,
                        no_weight_interpolation=True, all_markers_per_entry=True,
                        chronological_order='even-odd-to-power128', common_state_potential=False)
        try:
            with ctx.workprec(bank.precision):
                count = counts.count_diagnostics(cell, budget, cache=self.cache)
                identity['count'] = count
                result = _cell_bound(self.table, cell, self.prefix, count['support_count_upper'], budget)
                budget.check()
            identity.update(result, status='complete', complete_cell=True, directed=True,
                            meets_cell_budget=result['zero'] or F(result['log2_upper_rational']) <= -failure_bits)
        except (TimeoutError, ArithmeticError) as error:
            identity.update(failure=dict(error_type=type(error).__name__, message=str(error)),
                            meets_cell_budget=False)
        self.check_pins()
        return identity
