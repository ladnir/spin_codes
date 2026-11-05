"""Directed fixed-deficit covers of an explicit integer-weight band.

Every original integer weight is normalized through the frozen RegionalCache.
The weight-dependent outer factor is absorbed into all 256 regional factors
before maxima are taken. Raw midpoint rectangle sums count supports across
the band. A completed partition is required for any band probability bound.
"""

import ast
from collections import deque
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_gamma_all_pair_certificate as provider
import hd_field128_gamma_midpoint_boxes as boxes
import hd_field128_gamma_profile_certificate as profiles
import hd_field128_lift_alias_cover as restored


base, tail = provider.base, provider.tail
DIRECTORY = Path(__file__).resolve().parent
NOTE = '../HD_FIELD128_MIDPOINT_WEIGHT_BANDS.md'


def source_pins():
    pending, pins = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (DIRECTORY/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (DIRECTORY/target).is_file() and target not in pins:
                    pending.append(target)
    pins.update(boxes.source_pins())
    pins[NOTE] = sha256((DIRECTORY/NOTE).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    # No AST walks or closure discovery in the rectangle loop.
    if any(sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest
           for name, digest in pins.items()):
        raise RuntimeError('A pinned weight-band source changed')


def checked_weights(pairs, first, last):
    if (any(type(value) is not int for value in (pairs, first, last))
            or pairs < 1 or not 1 <= first <= last < 2*pairs or last-first+1 > 64):
        raise ValueError('Expected 1..64 consecutive weights in 1..2*pairs-1')
    return tuple(range(first, last+1))


def deficit_limit(pairs, weight):
    u, eps = divmod(weight, 2)
    return min(u, pairs-u-eps)


@dataclass(frozen=True)
class Rectangle:
    parity: int
    ranges: tuple

    def __post_init__(self):
        if (type(self.parity) is not int or self.parity not in (0, 1)
                or not isinstance(self.ranges, tuple) or len(self.ranges) != 3
                or any(not isinstance(row, tuple) or len(row) != 2
                       or any(type(x) is not int for x in row) or not 0 <= row[0] <= row[1]
                       for row in self.ranges)):
            raise ValueError('Expected parity and three nonnegative inclusive integer intervals U,D,E')

    def descriptor(self):
        return dict(parity=self.parity, u_deficit_ranges=[list(row) for row in self.ranges])

    def split(self, axis, cut):
        if (type(axis) is not int or axis not in (0, 1, 2) or type(cut) is not int
                or not self.ranges[axis][0] <= cut < self.ranges[axis][1]):
            raise ValueError('A split must leave two nonempty integer intervals')
        left, right = list(self.ranges), list(self.ranges)
        left[axis] = (self.ranges[axis][0], cut)
        right[axis] = (cut+1, self.ranges[axis][1])
        return Rectangle(self.parity, tuple(left)), Rectangle(self.parity, tuple(right))


def _maximum(left, right):
    return (max(left[0], right[0]), max(left[1], right[1]),
            max(left[2], right[2]), max(left[3], right[3]))


class RangeMax:
    """Immutable fixed-four-entry sparse table with O(1) range maxima."""

    def __init__(self, rows, budget):
        levels, half, size = [tuple(rows)], 1, len(rows)
        if not size:
            raise ValueError('Empty maximum table')
        while 2*half <= size:
            previous, following = levels[-1], []
            for start in range(size-2*half+1):
                if not start % 128:
                    budget.check()
                following.append(_maximum(previous[start], previous[start+half]))
            levels.append(tuple(following))
            half *= 2
        self.levels, self.last = tuple(levels), size-1

    def query(self, low, high):
        if (type(low) is not int or type(high) is not int
                or not 0 <= low <= high <= self.last):
            raise ValueError('Range outside prepared deficit domain')
        degree = (high-low+1).bit_length()-1
        rows, span = self.levels[degree], 1 << degree
        return _maximum(rows[low], rows[high-span+1])


class BandTables:
    """Generic arithmetic core, not by itself a production certificate.

    regional(weight,budget) must return all four-entry matrices at that
    original weight. Count prefixes use O(B*r) integers; O(B*r log r)
    directed matrices are retained. No profile-pair table is built.
    """

    def __init__(self, pairs, first, last, regional, constant, slope, *,
                 precision=256, repetitions=256, budget=None, progress=None):
        self.weights = checked_weights(pairs, first, last)
        if (type(precision) is not int or not 128 <= precision <= 4096
                or type(repetitions) is not int or repetitions < 2 or repetitions % 2):
            raise ValueError('Invalid precision or positive even region count')
        budget = budget or base.Budget()
        self.pairs, self.precision, self.repetitions = pairs, precision, repetitions
        count_tables, limits = {0: [], 1: []}, {}
        self.regional_histograms = 0
        with ctx.workprec(precision):
            if not constant.is_finite() or not slope.is_finite():
                raise ArithmeticError('Nonfinite outer prefix')
            self.constant = constant.upper()
            for weight in self.weights:
                budget.check()
                regions = regional(weight, budget)
                full = tail.full_counts(pairs, weight)
                if set(regions) != set(full):
                    raise ValueError('Every original-weight normalized histogram is required')
                u, eps = divmod(weight, 2)
                limit = deficit_limit(pairs, weight)
                target = [(arb(0),)*4 for unused in range(limit+1)]
                limits[eps] = max(limits.get(eps, 0), limit)
                factor = base.upper(arb(2)**(slope*weight/repetitions))
                for index, f in enumerate(full):
                    if not index % 128:
                        budget.check()
                    original = regions[f]
                    if len(original) != 4:
                        raise ValueError('Every regional matrix requires four entries')
                    row = tuple(base.upper(base.upper(value)*factor) for value in original)
                    d = u-f
                    target[d] = row
                # The raw midpoint rectangle does not use the costly near-face
                # deletion prefixes. A cutoff-zero count object suffices.
                counts = boxes.midpoint.MidpointCounts(pairs, weight, max_singletons=0, budget=budget)
                table = boxes.MidpointBoxes(pairs, weight, counts=counts, budget=budget)
                count_tables[eps].append((weight, u, limit, table, RangeMax(target, budget)))
                self.regional_histograms += len(full)
                if progress:
                    progress(dict(stage='prepared-original-weight', active_groups=weight,
                                  normalized_histograms=len(full)))
        self.count_tables = {eps: tuple(rows) for eps, rows in count_tables.items() if rows}
        self.limits = limits
        self.roots = tuple(Rectangle(eps, ((rows[0][1], rows[-1][1]),
                                           (0, limits[eps]), (0, limits[eps])))
                           for eps, rows in sorted(self.count_tables.items()))
        self.root_by_parity = {cell.parity: cell for cell in self.roots}
        self.prepared_weights = tuple(weight for weight in self.weights)
        budget.check()

    def validate_cell(self, cell):
        if type(cell) is not Rectangle or cell.parity not in self.root_by_parity:
            raise ValueError('Cell outside prepared parity/deficit domain')
        root = self.root_by_parity[cell.parity]
        if any(lo < outer[0] or hi > outer[1] for (lo, hi), outer in zip(cell.ranges, root.ranges)):
            raise ValueError('Cell outside prepared parity/deficit domain')

    def count(self, cell, budget):
        self.validate_cell(cell)
        (ul, uh), (dl, dh), (el, eh) = cell.ranges
        total = 0
        for weight, u, limit, table, maxima in self.count_tables[cell.parity]:
            budget.check()
            if ul <= u <= uh and dl <= limit and el <= limit:
                # rectangle() clips original full counts, so invalid deficits
                # at this particular weight contribute exactly zero.
                total += table.rectangle(u-dh, u-dl, u-eh, u-el)
        return total

    def matrices(self, cell, budget):
        self.validate_cell(cell)
        (ul, uh), (dl, dh), (el, eh) = cell.ranges
        left = right = (arb(0),)*4
        for weight, u, limit, table, maxima in self.count_tables[cell.parity]:
            budget.check()
            if ul <= u <= uh and dl <= limit and el <= limit:
                left = _maximum(left, maxima.query(dl, min(dh, limit)))
                right = _maximum(right, maxima.query(el, min(eh, limit)))
        return left, right

    def receipt(self, cell, failure_bits, budget):
        if type(failure_bits) is not int or not 0 <= failure_bits <= 4096:
            raise ValueError('Invalid cell failure budget')
        self.validate_cell(cell)
        with ctx.workprec(self.precision):
            count = self.count(cell, budget)
            transfer = arb(0)
            if count:
                left, right = self.matrices(cell, budget)
                if self.repetitions == 256:
                    transfer = profiles.pair_row_sum128(left, right)
                else:
                    matrix = base.matrix_power(base.matrix_product(left, right), self.repetitions//2)
                    transfer = base.upper(matrix[0]+matrix[1])
            zero = not count or transfer.is_zero()
            value = None if zero else self.constant+tail.log2arb(arb(count))+tail.log2arb(transfer)
            ceiling = None if zero else str(base.rational_ceiling(value))
            budget.check()
            return dict(cell=cell.descriptor(), status='complete', complete_cell=True,
                        count_upper=count, count_method='sum of raw midpoint rectangle caps over explicit original weights',
                        zero=zero, log2_upper_rational=ceiling,
                        log2_upper_enclosure=None if zero else str(value),
                        failure_bits=failure_bits,
                        meets_cell_budget=zero or F(ceiling) <= -failure_bits)


def _passes(receipt, cell, failure_bits, context_hash):
    if (receipt.get('cell') != cell.descriptor()
            or receipt.get('context_sha256') != context_hash):
        raise ValueError('Receipt cell/context mismatch')
    if receipt.get('complete_cell') is not True:
        return False
    if (receipt.get('status') != 'complete' or type(receipt.get('count_upper')) is not int
            or receipt['count_upper'] < 0 or receipt.get('failure_bits') != failure_bits):
        raise ValueError('Invalid completed cell receipt')
    if receipt.get('zero') is True:
        if receipt.get('log2_upper_rational') is not None:
            raise ValueError('Zero receipt has a finite logarithm')
        passing = True
    elif receipt.get('zero') is False and type(receipt.get('log2_upper_rational')) in (str, int):
        passing = F(receipt['log2_upper_rational']) <= -failure_bits
    else:
        raise ValueError('Missing explicit zero/logarithmic cell bound')
    if receipt.get('meets_cell_budget') is not passing:
        raise ValueError('Cell budget flag disagrees with ceiling')
    return passing


def audit_tree(nodes, roots, failure_bits, context_hash):
    by_id = {row['id']: row for row in nodes}
    if len(by_id) != len(nodes):
        raise ValueError('Duplicate tree node')
    accepted, pending, seen = [], [], set()

    def visit(identity, expected, parent):
        if identity in seen or identity not in by_id:
            raise ValueError('Missing, cyclic, or shared tree node')
        seen.add(identity)
        row = by_id[identity]
        if row['cell'] != expected.descriptor() or row['parent'] != parent:
            raise ValueError('Tree geometry/parent differs from exact split')
        if row['status'] == 'split':
            children = row['children']
            if len(children) != 2:
                raise ValueError('Both children must be retained')
            cells = expected.split(row['split']['axis'], row['split']['integer_cut'])
            for child, cell in zip(children, cells):
                visit(child, cell, identity)
        elif row['status'] == 'accepted':
            if not _passes(row['receipt'], expected, failure_bits, context_hash):
                raise ValueError('Accepted leaf lacks a passing complete receipt')
            accepted.append(identity)
        elif row['status'] in ('pending', 'unresolved'):
            pending.append(identity)
        else:
            raise ValueError('Unknown tree status')

    for index, root in enumerate(roots):
        visit(index, root, None)
    if seen != set(by_id):
        raise ValueError('Unreachable tree node')
    return accepted, pending


def cover(evaluator, *, max_nodes=2048, max_seconds=180., cell_failure_bits=64,
          aggregate_bits=54, checkpoint=None, context_checkpoint=None):
    """Bounded integer rectangle tree; partial checkpoints have no band bound.

    The generic injected-evaluator interface exists only for tiny tests.
    Only the actual Evaluator class produces a directed production result.
    Source guards run before and after the batch, including exceptional exits.
    """
    if (type(max_nodes) is not int or not 1 <= max_nodes <= 1000000
            or any(type(x) is not int or not 0 <= x <= 4096 for x in (cell_failure_bits, aggregate_bits))
            or any(callback is not None and not callable(callback) for callback in (checkpoint, context_checkpoint))):
        raise ValueError('Invalid cover limits or callbacks')
    budget = base.Budget(max_seconds)
    if getattr(evaluator, '_active', False):
        raise RuntimeError('A weight-band cover is already active')
    evaluator.check_pins()
    context = evaluator.context()
    context_raw = json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    context_hash = sha256(context_raw).hexdigest()
    roots = evaluator.tables.roots
    nodes = [dict(id=i, parent=None, cell=cell.descriptor(), status='pending') for i, cell in enumerate(roots)]
    cells, queue, blocked = list(roots), deque(range(len(roots))), []
    evaluations, sequence, stop_reason = 0, 0, None
    evaluator._active = True

    def publish(changed, action):
        nonlocal sequence
        if checkpoint:
            checkpoint(dict(sequence=sequence, action=action, context_sha256=context_hash,
                            changed_nodes=deepcopy(changed),
                            pending_frontier=list(queue)+blocked, evaluated_nodes=evaluations,
                            complete_requested_band=False, covers_all_supports_in_band=False,
                            directed=False, batch_source_verified_at_exit=False,
                            log2_upper_rational=None, whole_code_certificate=False))
        sequence += 1

    try:
        if context_checkpoint:
            context_checkpoint(dict(context=deepcopy(context), context_sha256=context_hash,
                                    directed=False, complete_requested_band=False, whole_code_certificate=False))
        publish(nodes, 'initial-roots')
        with ctx.workprec(evaluator.tables.precision):
            while queue:
                if evaluations >= max_nodes:
                    stop_reason = 'max_nodes'
                    break
                try:
                    budget.check()
                except TimeoutError:
                    stop_reason = 'max_seconds'
                    break
                index = queue.popleft()
                cell, node = cells[index], nodes[index]
                try:
                    receipt = evaluator.tables.receipt(cell, cell_failure_bits, budget)
                    receipt['context_sha256'] = context_hash
                except (TimeoutError, ArithmeticError) as error:
                    receipt = dict(cell=cell.descriptor(), context_sha256=context_hash, status='pending',
                                   complete_cell=False, log2_upper_rational=None,
                                   failure=dict(error_type=type(error).__name__, message=str(error)))
                evaluations += 1
                if _passes(receipt, cell, cell_failure_bits, context_hash):
                    nodes[index] = dict(node, status='accepted', receipt=receipt)
                    publish([nodes[index]], 'accepted-leaf')
                elif receipt['complete_cell'] is not True:
                    nodes[index] = dict(node, status='unresolved', receipt=receipt)
                    blocked.append(index)
                    publish([nodes[index]], 'unresolved-leaf')
                else:
                    root = evaluator.tables.root_by_parity[cell.parity]
                    # Balance relative widths so U can tighten cross-weight
                    # maxima before the much longer deficit axes are exhausted.
                    axis = max(range(3), key=lambda a: F(cell.ranges[a][1]-cell.ranges[a][0],
                                             max(1, root.ranges[a][1]-root.ranges[a][0])))
                    low, high = cell.ranges[axis]
                    if low == high:
                        nodes[index] = dict(node, status='unresolved', receipt=receipt)
                        blocked.append(index)
                        publish([nodes[index]], 'unsplittable-over-budget')
                    else:
                        cut = (low+high)//2
                        children = cell.split(axis, cut)
                        indices = (len(nodes), len(nodes)+1)
                        records = [dict(id=i, parent=index, cell=child.descriptor(), status='pending')
                                   for i, child in zip(indices, children)]
                        nodes[index] = dict(node, status='split', receipt=receipt,
                                            split=dict(axis=axis, integer_cut=cut), children=list(indices))
                        nodes.extend(records)
                        cells.extend(children)
                        queue.extend(indices)
                        publish([nodes[index]]+records, 'integer-split')
            accepted, pending = audit_tree(nodes, roots, cell_failure_bits, context_hash)
            complete = not pending
            result = dict(schema='gamma-midpoint-weight-band-cover-v1', context=context,
                          context_sha256=context_hash, status='complete' if complete else 'partial',
                          requested_weight_range=context['requested_weight_range'],
                          requested_weights=context['requested_weights'],
                          nodes=nodes, root_cells=[cell.descriptor() for cell in roots],
                          accepted_leaves=accepted, pending_frontier=pending, evaluated_nodes=evaluations,
                          tree_nodes=len(nodes), checkpoints=sequence,
                          complete_requested_band=complete, covers_all_supports_in_band=complete,
                          covers_all_matching_profiles=complete, every_integer_weight_prepared=True,
                          directed=complete and type(evaluator) is Evaluator,
                          numerical_backend='actual-all-pair-bank' if type(evaluator) is Evaluator else 'synthetic-test-only',
                          whole_code_certificate=False, no_weight_interpolation=True,
                          stop_reason=None if complete else stop_reason or 'unresolved-leaves',
                          zero=None, log2_upper_rational=None, log2_upper_enclosure=None,
                          meets_requested_band_target=False, requested_failure_bits=aggregate_bits,
                          cell_failure_bits=cell_failure_bits, max_nodes=max_nodes, max_seconds=max_seconds)
            if complete:
                total = sum((arb(2)**base.rational(F(nodes[i]['receipt']['log2_upper_rational']))
                             for i in accepted if not nodes[i]['receipt']['zero']), arb(0))
                zero = total.is_zero()
                value = None if zero else tail.log2arb(total)
                ceiling = None if zero else str(base.rational_ceiling(value))
                result.update(zero=zero, log2_upper_rational=ceiling,
                              log2_upper_enclosure=None if zero else str(value),
                              meets_requested_band_target=zero or F(ceiling) <= -aggregate_bits)
    finally:
        evaluator._active = False
        evaluator.check_pins()
    result.update(batch_source_verified_at_entry=True, batch_source_verified_at_exit=True,
                  elapsed_seconds=monotonic()-budget.started)
    return result


class Evaluator:
    """One actual bank and a fixed consecutive band of at most 64 weights."""

    def __init__(self, bank, first, last, *, max_seconds=180., witness_source=None, progress=None):
        if (type(bank) is not provider.AliasBank or bank.state != 30 or bank.q != 1 << 128
                or bank.distance != F(2, 5) or bank.outer_measure != 'invertible-toeplitz'):
            raise ValueError('The actual gamma30/d40 full-field Toeplitz all-pair AliasBank is required')
        checked_weights(base.PAIRS, first, last)
        self.bank, self.pins, self._active = bank, source_pins(), False
        self._witness_pin = None
        if witness_source is not None:
            if (not isinstance(witness_source, dict) or not isinstance(witness_source.get('path'), str)
                    or not isinstance(witness_source.get('sha256'), str)):
                raise ValueError('Witness provenance requires an explicit path and byte hash')
            path = Path(witness_source['path'])
            self._witness_pin = (path if path.is_absolute() else DIRECTORY/path, witness_source['sha256'])
            witness_source = deepcopy(witness_source)
        self.check_pins()
        budget = base.Budget(max_seconds)
        try:
            cache = restored.RegionalCache(bank)
            self.tables = BandTables(base.PAIRS, first, last, cache.regional,
                bank.constant, bank.slope, precision=bank.precision, budget=budget, progress=progress)
            context = dict(bank.metadata(), schema='gamma-midpoint-weight-band-context-v1',
                requested_weight_range=[first, last], requested_weights=list(self.tables.weights),
                parity_weights={str(eps): [r for r in self.tables.weights if r % 2 == eps]
                                for eps in sorted(self.tables.limits)},
                pairs=base.PAIRS, repetitions=256, paired_region_power=128,
                source_sha256=self.pins, bank_source_sha256=bank.pins, witness_source=witness_source,
                regional_histograms_prepared=self.tables.regional_histograms,
                normalized_every_original_weight=True,
                scaling='upward 2^(bank.slope*r/256) times each normalized regional entry; bank.constant once per event',
                max_domain='exact per-cell U range; both deficit intervals marginally nonempty; componentwise maxima of normalized scaled matrices',
                count_method='exact rectangle sums of raw midpoint caps, summed once over explicit weights',
                count_pointwise_minimum_used=False, count_only_complement_above_pairs=True,
                no_profile_pair_table=True, no_weight_interpolation=True,
                root_cells=[cell.descriptor() for cell in self.tables.roots],
                complete_requested_band=False, directed=False, whole_code_certificate=False,
                preparation_seconds=monotonic()-budget.started)
            self._context_raw = json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False)
        finally:
            self.check_pins()

    def check_pins(self):
        check_source_pins(self.pins)
        self.bank.check_pins()
        if self._witness_pin is not None:
            path, digest = self._witness_pin
            if sha256(path.read_bytes()).hexdigest() != digest:
                raise RuntimeError('The pinned witness source changed')

    def context(self):
        return json.loads(self._context_raw)

    def run_cover(self, **options):
        return cover(self, **options)
