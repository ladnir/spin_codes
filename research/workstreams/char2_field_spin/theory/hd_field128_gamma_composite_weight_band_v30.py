"""Fresh multi-witness band bounds over one exact shared integer partition.

Each frozen BandTables receipt already includes its witness's constant and
original-weight slope scaling. Minimize those COMPLETE cell bounds, never
their matrices or prefactors. The frozen v25 generic partition result retains
its original synthetic/directed flags verbatim. A separately identified
production envelope is issued only for this module's genuine fresh-bank
CompositeEvaluator. Existing single-bank auditors must not admit this schema.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import gc
from hashlib import sha256
import importlib
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma_profile_min_certificate_v29 as selected

HERE, NAME = Path(__file__).resolve().parent, Path(__file__).name
SELECTED_SHA = 'b79bd9edbeaab13f6adebb39cc2d73c16a784d453b54e96c73208c3f7d45d5e3'
CONTROLLER_NAME = 'hd_field128_gamma_portable_band_runner_v25'
CONTROLLER_SHA = '5d22c83e3419b1dc2a484b09c9cb2ada17971bd48195987140c9158942d7ef1c'
SCHEMA = 'gamma-composite-weight-band-production-v30'
MAX_NODES, COMPUTE_SECONDS, MAX_OUTPUT_BYTES = 4096, 540., 32 << 20


def require(test, message):
    if not test:
        raise ValueError(message)


def helpers():
    for name, expected in ((selected.NAME, SELECTED_SHA), (CONTROLLER_NAME+'.py', CONTROLLER_SHA)):
        require(sha256((HERE/name).read_bytes()).hexdigest() == expected, 'Frozen helper changed: '+name)
    selected.helpers()
    return importlib.import_module(CONTROLLER_NAME)


def source_pins(controller):
    pins = selected.source_pins(selected.helpers())
    for name, digest in controller.source_pins().items():
        require(name not in pins or pins[name] == digest, 'Conflicting source closure')
        pins[name] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def validate_request(first, last, max_nodes, cell_bits, aggregate_bits):
    require(type(first) is int and type(last) is int and 1 <= first <= last <= 8191
            and last-first+1 <= 64, 'Explicit 1..64 original weights in 1..8191 required')
    require(type(max_nodes) is int and 1 <= max_nodes <= MAX_NODES, 'Node limit must be 1..4096')
    require(type(cell_bits) is int and type(aggregate_bits) is int
            and 40 <= aggregate_bits <= cell_bits <= 4096, 'Expected 40 <= aggregate <= cell bits')


def load_proposals(specifications):
    """Read literal calibration witnesses; no saved bound or ancestry is imported."""
    require(2 <= len(specifications) <= 4, 'Expected 2..4 explicit pinned proposals')
    result, seen = [], set()
    for name, digest in specifications:
        name = selected.literals.proposal_name(name)
        require(name not in seen, 'Duplicate proposal path')
        require(isinstance(digest, str) and len(digest) == 64 and all(c in '0123456789abcdef' for c in digest),
                'Expected lowercase SHA-256')
        path = HERE/name
        require(path.resolve().is_relative_to(HERE.resolve()) and path.stat().st_size <= 1 << 20,
                'Escaping or oversized proposal')
        raw = path.read_bytes()
        require(sha256(raw).hexdigest() == digest, 'Proposal byte pin mismatch')
        saved = json.loads(raw)
        # This initial adapter accepts the reviewed lower-half calibration
        # format, but freshly normalizes every requested weight, including r>P.
        calibration_weight = selected.literals.validate_target(saved['active_groups'])
        witness, key = selected.literals.project_literals(saved, calibration_weight)
        result.append(dict(source=dict(path=name, sha256=digest), witness=witness,
            source_marker_key=key, calibration_weight=calibration_weight,
            source_role='Exact literals only; every requested original weight is normalized afresh'))
        seen.add(name)
    return result


class CompositeTables:
    """Generic tiny-testable arithmetic adapter; not a production credential."""
    def __init__(self, tables, bands):
        require(2 <= len(tables) <= 4 and all(type(t) is bands.BandTables for t in tables),
                'Expected 2..4 genuine frozen BandTables')
        first = tables[0]
        for table in tables[1:]:
            require(all(getattr(table, key) == getattr(first, key)
                for key in ('weights', 'pairs', 'precision', 'repetitions', 'roots', 'limits')),
                'Witnesses must bound identical original-weight/profile domains')
        self.components, self.bands = tuple(tables), bands
        self.weights, self.pairs = first.weights, first.pairs
        self.precision, self.repetitions = first.precision, first.repetitions
        self.roots, self.root_by_parity, self.limits = first.roots, first.root_by_parity, first.limits

    def receipt(self, cell, failure_bits, budget):
        receipts = []
        for table in self.components:
            budget.check()
            row = table.receipt(cell, failure_bits, budget)
            require(row['cell'] == cell.descriptor() and row['complete_cell'] is True
                and row['status'] == 'complete' and row['failure_bits'] == failure_bits,
                'Only complete bounds of the same cell may be minimized')
            require(type(row['count_upper']) is int and row['count_upper'] >= 0,
                    'Invalid integer cell count')
            receipts.append(row)
        require(all(row['count_upper'] == receipts[0]['count_upper']
                    and row['count_method'] == receipts[0]['count_method'] for row in receipts),
                'Component count domains differ')
        winner = min(range(len(receipts)), key=lambda i:
                     (not receipts[i]['zero'], F(0) if receipts[i]['zero'] else F(receipts[i]['log2_upper_rational'])))
        chosen = receipts[winner]
        result = dict(chosen, component_winner=winner,
            component_bounds=[dict(witness_index=i, zero=row['zero'],
                log2_upper_rational=row['log2_upper_rational'],
                log2_upper_enclosure=row['log2_upper_enclosure']) for i, row in enumerate(receipts)],
            minimum_stage='complete frozen cell bounds, including each own outer constant and original-r scaling',
            entrywise_cross_witness_minimum=False)
        result['meets_cell_budget'] = chosen['zero'] or F(chosen['log2_upper_rational']) <= -failure_bits
        return result


def prepare_component(proposal, first, last, controller, budget):
    """Build/check a real bank, detach frozen tables, then reclaim bank cycles."""
    bank = evaluator = None
    try:
        witness = proposal['witness']
        bank = controller.cert.AliasBank(F(witness['z']), F(witness['alpha']),
            [tuple(map(F, row)) for row in witness['markers']], state=30, q=1 << 128,
            distance=F(2, 5), grid=128, precision=256, bits=80,
            outer_measure='invertible-toeplitz', max_seconds=selected.remaining(budget),
            max_table_bytes=64 << 20)
        require(type(bank) is controller.cert.AliasBank, 'Actual frozen AliasBank required')
        evaluator = controller.bands.Evaluator(bank, first, last,
            max_seconds=selected.remaining(budget), witness_source=proposal['source'])
        evaluator.check_pins()
        tables, context = evaluator.tables, evaluator.context()
        require(type(tables) is controller.bands.BandTables
                and tables.weights == tuple(range(first, last+1)), 'Original-weight preparation incomplete')
        metadata = dict(bank.metadata(), exact_input_witness=deepcopy(witness),
            proposal_source=deepcopy(proposal['source']), calibration_weight=proposal['calibration_weight'],
            requested_weight_range=[first, last], fresh_bank_rebuilt=True,
            literal_source_role='calibration literals only; no calibration bound inherited')
        budget.check()
        return tables, context, metadata
    finally:
        if evaluator is not None:
            evaluator.check_pins()
        elif bank is not None:
            bank.check_pins()
        evaluator = bank = None
        gc.collect()


class CompositeEvaluator:
    """Production-only constructor: every component comes from a fresh actual bank."""
    def __init__(self, proposals, first, last, controller, budget, publish_component):
        validate_request(first, last, 1, 54, 48)
        require(2 <= len(proposals) <= 4, 'Expected 2..4 literal proposals')
        self.proposals, self.pins, self._active = deepcopy(proposals), source_pins(controller), False
        self._prepared = False
        components, lineage = [], []
        self.check_pins()
        try:
            for index, proposal in enumerate(self.proposals):
                tables, context, metadata = prepare_component(proposal, first, last, controller, budget)
                components.append(tables)
                self.check_pins()
                lineage.append(publish_component(index, context, metadata))
            self.tables = CompositeTables(components, controller.bands)
            context = dict(schema='gamma-composite-weight-band-context-v30',
                construction=selected.literals.CONSTRUCTION, field_order=1 << 128,
                state_dimension=30, target_distance='2/5', requested_weight_range=[first, last],
                requested_weights=list(self.tables.weights), proposals=self.proposals,
                components=lineage, source_sha256=self.pins, pairs=4096,
                precision=256, repetitions=256, paired_region_power=128,
                common_integer_count_domain=True, original_weights_normalized_individually=True,
                component_constants_and_slopes_kept_separate=True, entrywise_cross_witness_minimum=False,
                historical_origin_verified=False, old_numerical_bounds_reused=False,
                root_cells=[cell.descriptor() for cell in self.tables.roots],
                complete_requested_band=False, directed=False, whole_code_certificate=False)
            self._context_raw = json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False)
            self._prepared = True
        finally:
            self.check_pins()

    def context(self):
        require(self._prepared, 'Composite preparation is incomplete')
        return json.loads(self._context_raw)

    def check_pins(self):
        selected.check_pins(self.pins, self.proposals)


def production_envelope(evaluator, partition, controller):
    """New identity, never relabel the generic controller's nested result."""
    require(type(evaluator) is CompositeEvaluator and evaluator._prepared
        and type(evaluator.tables) is CompositeTables, 'Genuine prepared CompositeEvaluator required')
    evaluator.check_pins()
    context = evaluator.context()
    digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    require(partition['context'] == context and partition['context_sha256'] == digest,
            'Generic partition context changed')
    require(partition['directed'] is False and partition['numerical_backend'] == 'synthetic-test-only'
        and partition['whole_code_certificate'] is False, 'Generic controller flags must remain untouched')
    require(partition['controller'] == dict(path=controller.NAME, sha256=CONTROLLER_SHA), 'Wrong partition controller')
    accepted, pending = controller.bands.audit_tree(partition['nodes'], evaluator.tables.roots,
        partition['cell_failure_bits'], digest)
    require(accepted == partition['accepted_leaves'] and pending == partition['pending_frontier'],
            'Partition/frontier changed')
    complete = not pending
    require(partition['complete_requested_band'] is complete and partition['covers_all_matching_profiles'] is complete,
            'Generic completion contradicts exact tree')
    return dict(schema=SCHEMA, construction=selected.literals.CONSTRUCTION, field_order=1 << 128,
        state_dimension=30, target_distance='2/5', requested_weight_range=context['requested_weight_range'],
        requested_weights=context['requested_weights'], context_sha256=digest,
        composite_source_sha256=evaluator.pins, status='complete' if complete else 'partial',
        complete_requested_band=complete, covers_all_matching_profiles=complete, directed=complete,
        numerical_backend='fresh-actual-all-pair-banks-complete-cell-minimum',
        partition_result=partition, generic_partition_flags_preserved=True,
        zero=partition['zero'] if complete else None,
        log2_upper_rational=partition['log2_upper_rational'] if complete else None,
        log2_upper_enclosure=partition['log2_upper_enclosure'] if complete else None,
        requested_failure_bits=partition['requested_failure_bits'],
        meets_requested_band_target=complete and partition['meets_requested_band_target'],
        whole_code_certificate=False, historical_origin_verified=False,
        old_numerical_bounds_reused=False, no_weight_interpolation=True,
        scope='Separate fresh composite-band producer; nested generic engine is not the production credential')


def run(specifications, first, last, output_dir, *, max_nodes=4096, cell_bits=54, aggregate_bits=48):
    validate_request(first, last, max_nodes, cell_bits, aggregate_bits)
    with selected.remote_limits():
        controller = helpers()
        proposals = load_proposals(specifications)
        pins = source_pins(controller)
        check = lambda: selected.check_pins(pins, proposals)
        check()
        destination = HERE/controller.portable.relative_name(output_dir)
        controller.reserve_guard(HERE, MAX_OUTPUT_BYTES)
        destination.mkdir(exist_ok=False)
        artifacts = controller.OutputBudget(destination)
        budget, output_pins = controller.base.Budget(COMPUTE_SECONDS), []

        def put(name, value):
            path = destination/name
            artifacts.write_new(path, value)
            pin = dict(path=path.relative_to(HERE).as_posix(), sha256=sha256(path.read_bytes()).hexdigest())
            output_pins.append(pin)
            return pin

        def check_outputs():
            for pin in output_pins:
                require(sha256((HERE/pin['path']).read_bytes()).hexdigest() == pin['sha256'],
                        'Immutable composite output changed')

        def publish_component(index, context, metadata):
            return dict(index=index, context=put(f'component_{index:02d}_context.json', context),
                        bank_metadata=put(f'bank_{index:02d}.json', metadata))

        evaluator = None
        try:
            put('driver.json', dict(schema='gamma-composite-weight-band-driver-v30', source_sha256=pins,
                proposals=proposals, requested_weight_range=[first, last], maximum_nodes=max_nodes,
                cell_failure_bits=cell_bits, requested_failure_bits=aggregate_bits,
                hard_wall_seconds=600, cooperative_seconds=COMPUTE_SECONDS, memory_bytes=4 << 30,
                table_bytes=64 << 20, maximum_output_bytes=MAX_OUTPUT_BYTES,
                minimum_free_bytes=20 << 30, maximum_simultaneous_banks=1,
                fresh_banks_required=True, directed=False, whole_code_certificate=False))
            try:
                evaluator = CompositeEvaluator(proposals, first, last, controller, budget, publish_component)
                partition = controller.cover(evaluator, max_nodes=max_nodes,
                    max_seconds=selected.remaining(budget), cell_failure_bits=cell_bits, aggregate_bits=aggregate_bits,
                    checkpoint=lambda row: put(f"node_{row['sequence']:06d}.json", row),
                    context_checkpoint=lambda row: put('context.json', row))
                result = production_envelope(evaluator, partition, controller)
            except (TimeoutError, MemoryError, ArithmeticError) as error:
                result = dict(schema=SCHEMA, requested_weight_range=[first, last],
                    requested_weights=list(range(first, last+1)), status='incomplete',
                    complete_requested_band=False, covers_all_matching_profiles=False,
                    directed=False, whole_code_certificate=False, log2_upper_rational=None,
                    meets_requested_band_target=False,
                    failure=dict(error_type=type(error).__name__, message=str(error)),
                    composite_source_sha256=pins)
            check(); check_outputs()
            result.update(total_compute_seconds=monotonic()-budget.started,
                          bytes_before_summary=artifacts.used_bytes,
                          runner=dict(path=NAME, sha256=pins[NAME]))
            summary = put('summary.json', result)
            return dict(status=result['status'], requested_weight_range=[first, last], summary=summary,
                complete_requested_band=result['complete_requested_band'],
                log2_upper_rational=result['log2_upper_rational'],
                meets_requested_band_target=result['meets_requested_band_target'],
                output_bytes=artifacts.used_bytes, whole_code_certificate=False)
        finally:
            if evaluator is not None:
                evaluator.check_pins()
            check(); check_outputs()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', nargs=2, action='append', required=True, metavar=('RELATIVE_JSON', 'SHA256'))
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--max-nodes', type=int, default=4096)
    parser.add_argument('--cell-bits', type=int, default=54)
    parser.add_argument('--aggregate-bits', type=int, default=48)
    args = parser.parse_args()
    result = run(args.proposal, args.first, args.last, args.output_dir,
        max_nodes=args.max_nodes, cell_bits=args.cell_bits, aggregate_bits=args.aggregate_bits)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 0 if result['complete_requested_band'] and result['meets_requested_band_target'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
