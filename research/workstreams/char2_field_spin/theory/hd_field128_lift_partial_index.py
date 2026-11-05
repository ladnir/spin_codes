"""Index explicitly supplied, compatible literal-GF128 probability receipts.

For weight r, the indexed event is the existence of a forbidden output from
any nonzero input having r active outer groups. Each retained row must cover
all supports and all even/odd splits at that weight. A terminal receipt
covers the union of these events for 5854..8192; it enters the sum only once.

This tool checks provenance, metadata, and scalar enclosure consistency.
It does not replay regional coefficients or authenticate supplied receipts.
Missing weights are never interpolated. The output is always a partial
index, even if its coverage and aggregate meet the whole-code thresholds.
"""

import ast
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from glob import glob
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path

from flint import arb, ctx


DIRECTORY = Path(__file__).resolve().parent
CONSTRUCTION = 'toeplitz-fullfield-h4-alternating-swapped-pairs-k20-rate-half-v1'
FIELD_ORDER, GROUPS, PAIRS = 1 << 128, 8192, 4096
DIMENSION, OUTPUT_LENGTH = 1 << 20, 1 << 21
TERMINAL_RANGE = [5854, 8192]
MAX_INPUT_BYTES = 64 << 20
CERTIFICATE = 'hd_field128_lift_certificate.py'
SPARSE = 'hd_field128_lift_sparse_repair.py'
WIDE = 'hd_field128_lift_sparse_wide.py'
BANDS = 'hd_field128_lift_sparse_bands.py'
COVER = 'hd_field128_lift_alias_cover.py'
BATCH = 'hd_field128_lift_alias_batch.py'
TAIL = 'hd_field128_lift_tail_certificate.py'
ALIAS_METHODS = {
    'positive-residue-aliases-directed-interval-radix2-DFT':
        'hd_field128_lift_alias_certificate.py',
    'positive-residue-aliases-directed-native-ACB-DFT':
        'hd_field128_lift_acb_alias_certificate.py',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same(record, expected, label):
    require(isinstance(record, dict), label + ' must be an object')
    for key, value in expected.items():
        require(key in record and type(record[key]) is type(value) and record[key] == value,
                label + ': incompatible ' + key)


def fraction(value):
    require(isinstance(value, str) and 0 < len(value) <= 4096,
            'Expected a bounded exact rational string')
    return F(value)


def interval(value):
    require(isinstance(value, list) and len(value) == 2
            and all(type(x) is int for x in value)
            and 1 <= value[0] <= value[1] <= GROUPS, 'Invalid weight interval')
    return value


def ranges(weights):
    result = []
    for weight in sorted(set(weights)):
        if result and result[-1][1] + 1 == weight:
            result[-1][1] = weight
        else:
            result.append([weight, weight])
    return result


def range_weights(value):
    require(isinstance(value, list) and bool(value), 'Missing explicit weight ranges')
    weights = [weight for item in value for weight in range(interval(item)[0], item[1] + 1)]
    require(weights == sorted(set(weights)), 'Requested ranges overlap or are out of order')
    return weights


def source_closure(entrypoint):
    """Hash the recursive static local-import closure without importing it."""
    pending, result = [entrypoint], {}
    while pending:
        name = pending.pop()
        if name in result:
            continue
        require(Path(name).name == name and name.endswith('.py'), 'Invalid source name')
        path = DIRECTORY / name
        raw = path.read_bytes()
        result[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw, filename=str(path))):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0] + '.py'
                if (DIRECTORY / target).is_file() and target not in result:
                    pending.append(target)
    return dict(sorted(result.items()))


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError('Nonfinite JSON literal: ' + value)


def rational_arb(value):
    return arb(value.numerator) / value.denominator


def scalar_ceiling(row, precision):
    require(type(precision) is int and 128 <= precision <= 4096, 'Invalid scalar precision')
    bound = fraction(row.get('log2_upper_rational'))
    text = row.get('log2_upper_enclosure')
    require(isinstance(text, str) and 0 < len(text) <= 16384, 'Missing bounded scalar enclosure')
    with ctx.workprec(precision):
        value = arb(text)
        require(value.is_finite() and rational_arb(bound) >= value.upper(),
                'Scalar ceiling does not dominate its enclosure')
    return bound


def union_ceiling(bounds, precision=256):
    """Enclose log2(sum(2**b)); the input b values are exact rational ceilings."""
    if not bounds:
        return None
    with ctx.workprec(precision):
        total = sum((arb(2) ** rational_arb(bound) for bound in bounds), arb(0))
        value = total.log() / arb(2).log()
        require(value.is_finite(), 'Nonfinite indexed union')
        return F(int((value.upper() * 10**12).upper().ceil().fmpz()), 10**12)


@dataclass(frozen=True)
class Candidate:
    weight: int | None
    bound: F
    provenance: dict

    def record(self):
        return dict(active_groups=self.weight, log2_upper_rational=str(self.bound),
                    provenance=self.provenance)


class Validator:
    def __init__(self, distance):
        require(isinstance(distance, F) and 0 < distance < F(1, 2), 'Invalid exact target distance')
        self.distance = distance
        self.files = {}
        self.sources = {}
        self.closures = {}

    def common(self):
        return dict(construction=CONSTRUCTION, field_order=FIELD_ORDER, field_bits=128,
                    target_distance=str(self.distance), forbidden_weight=int(self.distance * OUTPUT_LENGTH),
                    outer_measure='invertible-toeplitz', directed=True, whole_code_certificate=False)

    def full(self):
        return dict(self.common(), source_dimension=DIMENSION, output_length=OUTPUT_LENGTH,
                    outer_groups=GROUPS, projective_divisor=str(FIELD_ORDER - 1), pair_swaps=True,
                    state_persists_across_regions=True, full_field_random_coefficients=True)

    def semantics(self, record, *, full=True):
        same(record, self.full() if full else self.common(), 'receipt semantics')
        # A format may inherit missing metadata, but contradictory metadata is never ignored.
        same(record, {key: value for key, value in self.full().items() if key in record},
             'present receipt semantics')

    def pins(self, record, key, producer):
        if producer not in self.closures:
            closure = source_closure(producer)
            for name, digest in closure.items():
                require(name not in self.sources or self.sources[name] == digest,
                        'Source changed between closure checks: ' + name)
                self.sources[name] = digest
            self.closures[producer] = closure
        same(record, {key: self.closures[producer]}, 'producer source pins')

    def load(self, path, *, expected=None, role='explicit-input'):
        path = Path(path).resolve()
        require(path.is_file() and path.stat().st_size <= MAX_INPUT_BYTES, 'Missing or oversized receipt: ' + str(path))
        raw = path.read_bytes()
        require(len(raw) <= MAX_INPUT_BYTES, 'Oversized receipt')
        digest = sha256(raw).hexdigest()
        if expected is not None:
            require(isinstance(expected, str) and digest == expected, 'Receipt hash mismatch: ' + str(path))
        prior = self.files.get(str(path))
        require(prior is None or prior['sha256'] == digest, 'Receipt changed while indexing: ' + str(path))
        pin = dict(path=str(path), sha256=digest, bytes=len(raw))
        self.files[str(path)] = dict(pin, roles=sorted(set((prior or {}).get('roles', []) + [role])))
        value = json.loads(raw.decode('utf-8'), object_pairs_hook=_unique_object, parse_constant=_reject_constant)
        require(isinstance(value, dict), 'Receipt must be a JSON object')
        return value, pin

    def dependency(self, reference, parent, role):
        require(isinstance(reference, dict) and isinstance(reference.get('path'), str), 'Missing dependency path')
        require(isinstance(reference.get('sha256'), str), 'Missing dependency hash')
        path = Path(reference['path'])
        if not path.is_absolute():
            path = Path(parent['path']).parent / path
        return self.load(path, expected=reference['sha256'], role=role)

    def row(self, row, parent, index, kind, precision):
        require(isinstance(row, dict), 'Weight row must be an object')
        weight = row.get('active_groups')
        require(type(weight) is int and 1 <= weight <= GROUPS, 'Invalid active-group weight')
        same(row, dict(covers_all_supports_at_weight=True, covers_all_even_odd_splits=True), 'weight coverage')
        same(row, {key: value for key, value in self.full().items() if key in row}, 'present row semantics')
        return Candidate(weight, scalar_ceiling(row, precision),
                         dict(parent, kind=kind, row_index=index))

    def alias_rows(self, receipt, pin):
        self.semantics(receipt)
        method = receipt.get('coefficient_method')
        require(method in ALIAS_METHODS, 'Unknown directed alias producer')
        self.pins(receipt, 'source_sha256', ALIAS_METHODS[method])
        self.pins(receipt, 'coverage_driver_sha256', COVER)
        require(receipt.get('support_method') in ('runs', 'marginal'), 'Unknown support summation method')
        witness = receipt.get('witness')
        require(isinstance(witness, dict), 'Missing alias witness')
        require(0 < fraction(witness.get('z')) < 1 and 0 < fraction(witness.get('alpha')) <= 1,
                'Invalid alias witness')
        precision, bits = witness.get('precision'), witness.get('packet_mantissa_bits')
        require(type(bits) is int and 32 <= bits <= 512, 'Invalid packet mantissa budget')
        rows = receipt.get('weight_receipts')
        require(isinstance(rows, list) and len(rows) <= GROUPS, 'Invalid alias rows')
        lo, hi = interval(receipt.get('requested_weight_range'))
        require(len(rows) <= hi - lo + 1, 'Too many alias rows')
        expected_weights = list(range(lo, lo + len(rows)))
        require([row.get('active_groups') for row in rows] == expected_weights, 'Alias rows are not the completed prefix')
        complete = len(rows) == hi - lo + 1
        same(receipt, dict(weights_evaluated=len(rows), every_integer_weight_evaluated=complete,
                           complete_requested_range=complete,
                           completed_weight_range=[lo, lo + len(rows) - 1] if rows else None,
                           pending_weight_ranges=[] if complete else [[lo + len(rows), hi]]), 'alias coverage counts')
        result = []
        for index, row in enumerate(rows):
            weight = row['active_groups']
            same(row, dict(regional_histograms=weight // 2 - max(0, weight - PAIRS) + 1), 'matching histogram coverage')
            require(fraction(row.get('state_potential_ratio')) > 0, 'Invalid state potential')
            result.append(self.row(row, pin, index, 'alias-cover', precision))
        return result

    def sparse_rows(self, receipt, pin):
        self.semantics(receipt)
        self.pins(receipt, 'source_sha256', SPARSE)
        self.pins(receipt, 'repair_source_sha256', CERTIFICATE)
        rows = receipt.get('selected_weight_receipts')
        require(isinstance(rows, list), 'Missing selected sparse rows')
        lo, hi = interval(receipt.get('requested_weight_range'))
        require(hi <= 64 and [row.get('active_groups') for row in rows] == list(range(lo, hi + 1)),
                'Sparse selected rows do not cover their declared interval')
        same(receipt, dict(selected_weight_count=len(rows), duplicate_weight_count=0,
                           complete_requested_range=True, every_weight_evaluated=True,
                           covers_all_supports_in_range=True, covers_all_even_odd_splits=True), 'sparse coverage')
        original, original_pin = self.dependency(receipt.get('input_provenance'), pin, 'sparse-original')
        require(original == receipt.get('input_receipt'), 'Embedded sparse input differs from its pinned file')
        self.semantics(original, full=False)
        self.pins(original, 'source_sha256', BANDS)
        completed = receipt.get('completed_repairs')
        require(isinstance(completed, list), 'Missing completed repair provenance')
        result = []
        for index, row in enumerate(rows):
            origin = row.get('origin')
            require(isinstance(origin, dict), 'Missing sparse row origin')
            clean = {key: value for key, value in row.items() if key != 'origin'}
            if origin.get('kind') == 'input-band-row':
                band_index, row_index = origin.get('band_index'), origin.get('row_index')
                bands = original.get('completed_bands')
                require(isinstance(bands, list) and type(band_index) is int and 0 <= band_index < len(bands),
                        'Invalid sparse band origin')
                band = bands[band_index]
                self.semantics(band)
                self.pins(band, 'source_sha256', BANDS)
                children = band.get('weight_receipts')
                require(isinstance(children, list) and type(row_index) is int and 0 <= row_index < len(children),
                        'Invalid sparse row origin')
                require(clean == children[row_index], 'Selected sparse row differs from its origin')
            elif origin.get('kind') == 'own-weight-repair':
                repair_index = origin.get('repair_index')
                require(type(repair_index) is int and 0 <= repair_index < len(completed), 'Invalid repair origin')
                require(clean == completed[repair_index], 'Selected sparse repair differs from its origin')
                self.semantics(clean)
                self.pins(clean, 'source_sha256', CERTIFICATE)
            else:
                raise ValueError('Unknown sparse row origin')
            same(row, dict(directed=True, whole_code_certificate=False), 'sparse row direction')
            witness = row.get('witness')
            require(isinstance(witness, dict), 'Missing sparse row witness')
            require(0 < fraction(witness.get('z')) < 1 and 0 < fraction(witness.get('alpha')) <= 1,
                    'Invalid sparse witness')
            same(witness, dict(paired=True), 'sparse witness')
            maximum = witness.get('coefficient_maximum')
            require(type(maximum) is int and row['active_groups'] <= maximum <= 64, 'Invalid sparse coefficient maximum')
            same(witness, dict(coefficient_marker=str(F(maximum, PAIRS))), 'sparse coefficient marker')
            candidate = self.row(row, pin, index, 'sparse-repaired', witness.get('precision'))
            result.append(Candidate(candidate.weight, candidate.bound,
                                    dict(candidate.provenance, row_origin=origin, original_input=original_pin)))
        return result

    def batch_rows(self, receipt, pin):
        self.semantics(receipt, full=False)
        same(receipt, dict(schema='field128-lift-alias-batch-v1', pair_swaps=True,
                           full_field_random_coefficients=True, old_numerical_ceilings_reused=False), 'batch semantics')
        self.pins(receipt, 'source_sha256', BATCH)
        lo, hi = interval(receipt.get('requested_weight_range'))
        rows, bands = receipt.get('weight_receipts'), receipt.get('completed_bands')
        require(isinstance(rows, list) and len(rows) <= GROUPS and isinstance(bands, list), 'Invalid batch rows/bands')
        children, expected = {}, []
        for band in bands:
            require(isinstance(band, dict) and isinstance(band.get('id'), str), 'Invalid completed band')
            child, child_pin = self.dependency(band.get('receipt'), pin, 'batch-band')
            child_rows = self.alias_rows(child, child_pin)
            child_lo, child_hi = interval(child.get('requested_weight_range'))
            require(lo <= child_lo <= child_hi <= hi, 'Batch child exceeds its requested interval')
            same(band, dict(weight_range=[child_lo, child_hi], complete_requested_range=child['complete_requested_range'],
                            completed_union_log2_upper_rational=child['completed_union_log2_upper_rational']), 'batch child metadata')
            for row, candidate in zip(child['weight_receipts'], child_rows):
                annotated = dict(row, band_id=band['id'], receipt_sha256=child_pin['sha256'])
                expected.append(annotated)
                key = (band['id'], child_pin['sha256'], candidate.weight)
                require(key not in children, 'Duplicate completed band row')
                children[key] = candidate
        require(rows == expected, 'Batch rows differ from the pinned child rows')
        covered = {row['active_groups'] for row in rows}
        missing = set(range(lo, hi + 1)) - covered
        same(receipt, dict(covered_weight_count=len(covered), duplicate_weight_count=len(rows) - len(covered),
                           pending_weight_ranges=ranges(missing), complete_requested_range=not missing), 'batch coverage counts')
        result = []
        for index, row in enumerate(rows):
            candidate = children[(row['band_id'], row['receipt_sha256'], row['active_groups'])]
            result.append(Candidate(candidate.weight, candidate.bound,
                                    dict(pin, kind='alias-batch', row_index=index, child=candidate.provenance)))
        return result

    def wide_band_rows(self, receipt, pin):
        self.semantics(receipt)
        self.pins(receipt, 'source_sha256', WIDE)
        same(receipt, dict(coefficient_method='positive-carry-free-Arb-truncated-bivariate-cutoff256',
                           covers_all_supports_at_completed_weights=True, covers_all_even_odd_splits=True),
             'wide sparse band')
        requested = range_weights(receipt.get('requested_weight_ranges'))
        require(requested[-1] <= 256, 'Wide sparse range exceeds supported cutoff')
        rows = receipt.get('weight_receipts')
        require(isinstance(rows, list) and len(rows) <= len(requested), 'Invalid wide sparse rows')
        require([row.get('active_groups') for row in rows] == requested[:len(rows)], 'Wide rows are not the completed prefix')
        missing = set(requested[len(rows):])
        same(receipt, dict(weights_evaluated=len(rows), every_requested_weight_evaluated=not missing,
                           complete_requested_ranges=not missing, pending_weight_ranges=ranges(missing)), 'wide band coverage')
        maximum = receipt.get('coefficient_maximum')
        require(type(maximum) is int and requested[-1] <= maximum <= 256, 'Invalid wide coefficient maximum')
        same(receipt, dict(coefficient_marker=str(F(maximum, PAIRS))), 'wide coefficient marker')
        witness = receipt.get('witness')
        require(isinstance(witness, dict), 'Missing wide witness')
        require(0 < fraction(witness.get('z')) < 1 and 0 < fraction(witness.get('alpha')) <= 1, 'Invalid wide witness')
        same(witness, dict(paired=True), 'wide witness')
        result = []
        for index, row in enumerate(rows):
            same(row, dict(directed=True, whole_code_certificate=False,
                           witness=dict(witness, coefficient_maximum=maximum, coefficient_marker=str(F(maximum, PAIRS)))),
                 'wide row witness')
            result.append(self.row(row, pin, index, 'sparse-wide-band', witness.get('precision')))
        return result

    def wide_rows(self, receipt, pin):
        self.semantics(receipt)
        self.pins(receipt, 'source_sha256', WIDE)
        plan, bands, rows = receipt.get('plan'), receipt.get('band_receipts'), receipt.get('weight_receipts')
        require(isinstance(plan, dict) and isinstance(bands, list) and isinstance(rows, list), 'Invalid wide sparse plan')
        same(plan, dict(old_numerical_ceilings_reused=False, directed=False, whole_code_certificate=False), 'wide proposal scope')
        requested = range_weights(receipt.get('requested_weight_ranges'))
        same(plan, dict(requested_weight_ranges=receipt['requested_weight_ranges']), 'wide plan ranges')
        require(requested[-1] <= 256, 'Wide plan exceeds supported cutoff')
        jobs = plan.get('jobs')
        require(isinstance(jobs, list) and len(bands) <= len(jobs), 'Invalid wide plan jobs')
        expected, result = [], []
        for index, band in enumerate(bands):
            job = jobs[index]
            require(isinstance(job, dict), 'Invalid wide job')
            same(job, dict(old_numerical_ceilings_reused=False), 'wide job proposal scope')
            same(band, dict(old_witness_job=index, requested_weight_ranges=job.get('requested_weight_ranges'),
                            coefficient_maximum=job.get('coefficient_maximum')), 'wide band job')
            same(band.get('witness'), dict(z=job.get('z'), alpha=job.get('alpha')), 'wide job witness')
            candidates = self.wide_band_rows(band, pin)
            for row, candidate in zip(band['weight_receipts'], candidates):
                expected.append(dict(row, job_index=index))
                result.append(Candidate(candidate.weight, candidate.bound,
                                        dict(candidate.provenance, kind='sparse-wide', band_index=index,
                                             row_index=len(result))))
        require(rows == expected, 'Wide plan rows differ from embedded band rows')
        covered = {row['active_groups'] for row in rows}
        require(covered <= set(requested), 'Wide plan row outside requested ranges')
        missing = set(requested) - covered
        same(receipt, dict(completed_weight_count=len(rows), duplicate_weight_count=len(rows) - len(covered),
                           complete_requested_ranges=not missing, pending_weight_ranges=ranges(missing)), 'wide plan coverage')
        return result

    def terminal(self, receipt, pin):
        self.semantics(receipt, full=False)
        same(receipt, dict(weight_range=TERMINAL_RANGE, source_dimension=DIMENSION, output_length=OUTPUT_LENGTH,
                           projective_divisor=str(FIELD_ORDER - 1), covers_all_supports_in_range=True,
                           covers_all_even_odd_splits=True, no_support_subset_union=True,
                           method='inverse-image-terminal-range'), 'terminal event')
        self.pins(receipt, 'source_sha256', TAIL)
        witness = receipt.get('witness')
        require(isinstance(witness, dict), 'Missing terminal witness')
        require(0 < fraction(witness.get('mark')) < 1, 'Invalid terminal mark')
        return Candidate(None, scalar_ceiling(receipt, witness.get('precision')),
                         dict(pin, kind='terminal-range', weight_range=TERMINAL_RANGE))

    def extract(self, receipt, pin):
        if receipt.get('method') == 'inverse-image-terminal-range':
            return [], [self.terminal(receipt, pin)]
        if receipt.get('schema') == 'field128-lift-alias-batch-v1':
            return self.batch_rows(receipt, pin), []
        if 'selected_weight_receipts' in receipt:
            return self.sparse_rows(receipt, pin), []
        if 'band_receipts' in receipt and 'requested_weight_ranges' in receipt:
            return self.wide_rows(receipt, pin), []
        if receipt.get('coefficient_method') == 'positive-carry-free-Arb-truncated-bivariate-cutoff256':
            return self.wide_band_rows(receipt, pin), []
        if 'coverage_driver_sha256' in receipt and 'weight_receipts' in receipt:
            return self.alias_rows(receipt, pin), []
        raise ValueError('Unsupported receipt; accepted formats: sparse-repaired, sparse-wide, alias-cover, alias-batch, terminal')

    def unchanged(self):
        for name, digest in self.sources.items():
            require(sha256((DIRECTORY / name).read_bytes()).hexdigest() == digest, 'Source changed during index: ' + name)
        for pin in self.files.values():
            require(sha256(Path(pin['path']).read_bytes()).hexdigest() == pin['sha256'], 'Receipt changed during index: ' + pin['path'])


def build_index(paths, *, distance=F(7, 20), target_bits=80, aggregate_bits=40, precision=256):
    require(type(target_bits) is int and 0 <= target_bits <= 4096, 'Invalid event target')
    require(type(aggregate_bits) is int and 0 <= aggregate_bits <= 4096, 'Invalid aggregate target')
    require(type(precision) is int and 128 <= precision <= 4096, 'Invalid aggregation precision')
    paths = sorted({str(Path(path).resolve()) for path in paths})
    require(bool(paths), 'Supply at least one explicit input or pattern')
    validator = Validator(distance)
    own_pins = source_closure(Path(__file__).name)
    validator.sources.update(own_pins)
    rows, terminals, inputs = [], [], []
    for path in paths:
        receipt, pin = validator.load(path)
        new_rows, new_terminals = validator.extract(receipt, pin)
        rows.extend(new_rows)
        terminals.extend(new_terminals)
        inputs.append(dict(pin, indexed_rows=len(new_rows), indexed_terminal_events=len(new_terminals)))
    key = lambda candidate: (candidate.bound, candidate.provenance['sha256'], candidate.provenance['path'],
                             candidate.provenance.get('row_index', -1))
    selected = {}
    counts = Counter(row.weight for row in rows)
    for row in rows:
        if row.weight not in selected or key(row) < key(selected[row.weight]):
            selected[row.weight] = row
    terminal = min(terminals, key=key) if terminals else None
    terminal_passes = terminal is not None and terminal.bound <= -target_bits
    terminal_weights = set(range(TERMINAL_RANGE[0], TERMINAL_RANGE[1] + 1)) if terminal_passes else set()
    passing = {weight for weight, row in selected.items() if row.bound <= -target_bits}
    # The terminal range is one union event, never 2339 separate probability terms.
    contributing_weights = sorted(passing - terminal_weights)
    events = [selected[weight].record() for weight in contributing_weights]
    if terminal_passes:
        events.append(dict(terminal.record(), weight_range=TERMINAL_RANGE))
    aggregate = union_ceiling([selected[weight].bound for weight in contributing_weights]
                              + ([terminal.bound] if terminal_passes else []), precision)
    covered = passing | terminal_weights
    universe = set(range(1, GROUPS + 1))
    known = set(selected) | (set(range(TERMINAL_RANGE[0], TERMINAL_RANGE[1] + 1)) if terminal else set())
    missing, unresolved = universe - known, universe - covered
    above_target = sorted(unresolved - missing)
    validator.unchanged()
    return dict(validator.full(), schema='field128-lift-partial-index-v1',
                status='partial receipt index; not a whole-code certificate',
                indexed_weight_domain=[1, GROUPS], individual_event_failure_bits=target_bits,
                aggregate_failure_bits=aggregate_bits, aggregation_precision=precision,
                passing_weight_ranges=ranges(covered), passing_weight_count=len(covered),
                unresolved_weight_ranges=ranges(unresolved), unresolved_weights=sorted(unresolved),
                missing_weight_ranges=ranges(missing), missing_weights=sorted(missing),
                known_above_target_weight_ranges=ranges(above_target), known_above_target_weights=above_target,
                selected_per_weight=[dict(selected[weight].record(), compatible_candidate_count=counts[weight],
                                          passes_individual_target=selected[weight].bound <= -target_bits,
                                          used_in_union=weight in contributing_weights)
                                     for weight in sorted(selected)],
                terminal_event=(dict(terminal.record(), weight_range=TERMINAL_RANGE,
                                     compatible_candidate_count=len(terminals), passes_individual_target=terminal_passes,
                                     used_in_union=terminal_passes) if terminal else None),
                union_events=events, union_event_count=len(events), terminal_union_event_count=int(terminal_passes),
                selected_union_log2_upper_rational=str(aggregate) if aggregate is not None else None,
                all_weights_covered=not unresolved,
                selected_union_meets_aggregate_target=aggregate is not None and aggregate <= -aggregate_bits,
                whole_code_threshold_conditions_met=(not unresolved and aggregate is not None and aggregate <= -aggregate_bits),
                explicit_inputs=inputs, receipt_file_pins=sorted(validator.files.values(), key=lambda pin: pin['path']),
                source_sha256=own_pins, checked_producer_closures=validator.closures,
                checked_source_sha256=dict(sorted(validator.sources.items())),
                runtime=dict(python_flint=version('python-flint')),
                source_pins_unchanged_at_exit=True, input_files_unchanged_at_exit=True,
                selection_rule='minimum exact Fraction ceiling for each compatible event; deterministic hash/path ties',
                terminal_rule='best compatible terminal receipt; one term replaces all per-weight terms in 5854..8192 when it passes',
                input_validation='construction, directed metadata, source closures, pinned row origins, scalar ceilings; no coefficient replay',
                provenance_scope='hashes record file agreement and provenance, not receipt authentication',
                dependency_scope='explicit inputs plus their named sparse-original and batch-band dependencies; no ensemble discovery',
                scope='union of listed passing events only; no claim for unresolved weights; whole_code_certificate remains false')


def explicit_paths(inputs, patterns):
    result = list(inputs)
    for pattern in patterns:
        matches = sorted(glob(pattern, recursive=False))
        require(bool(matches), 'Explicit pattern matched no files: ' + pattern)
        result.extend(matches)
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', action='append', default=[], help='Explicit receipt path; repeat as needed')
    parser.add_argument('--pattern', action='append', default=[], help='Explicit nonrecursive input glob; no default patterns')
    parser.add_argument('--distance', default='7/20', help='Exact target distance; default 7/20')
    parser.add_argument('--target-bits', type=int, default=80)
    parser.add_argument('--aggregate-bits', type=int, default=40)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', required=True, type=Path, help='New JSON output path; exclusive create')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists; choose a new path')
    result = build_index(explicit_paths(args.input, args.pattern), distance=F(args.distance),
                         target_bits=args.target_bits, aggregate_bits=args.aggregate_bits, precision=args.precision)
    result['input_selection'] = dict(explicit_paths=args.input, explicit_patterns=args.pattern)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('status', 'passing_weight_ranges', 'unresolved_weight_ranges',
                                                  'union_event_count', 'selected_union_log2_upper_rational',
                                                  'all_weights_covered', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
