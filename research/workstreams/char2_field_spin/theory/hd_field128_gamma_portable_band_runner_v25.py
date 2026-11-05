"""Compact absolute-width band controller; frozen numerical arithmetic.

The evaluator, cell bounds, integer splits, tree audit and positive union law
are unchanged. Only partition scheduling and checkpoint storage differ from
v24: ties prefer D, then E, then U; checkpoints omit pending_frontier, while
the final summary retains the complete tree and frontier. This is not a new
construction. The legacy cover schema is the numerical receipt wire format;
the separate controller identity and source closure identify this producer.

Production is a Linux-only bounded run for assigned weights 368..383. It uses
a 4 GiB address-space ceiling, a shared cooperative 600-second compute budget,
32 MiB of JSON output, and a post-write 20 GiB disk reserve. No saved numerical
bound is consumed. Only the portable seed's exact positive literals are used.
"""
from collections import deque
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import argparse
import json
from math import isfinite
from pathlib import Path
import shutil
import sys
from time import monotonic

from flint import arb, ctx

import hd_field128_gamma_midpoint_weight_bands as bands
import hd_field128_gamma_portable_band_runner_v24 as portable


HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FROZEN_BANDS = '4ce8c0938379f267cc75097205fe2ece5616474e5db418efef4f153e02e068a8'
FROZEN_PORTABLE = 'e69c97a06032bbe06aade7d9ed707f9353b633c464e70d740d8647ddea58ed1b'
MAX_NODES, MAX_SECONDS = 4096, 600.
MAX_OUTPUT_BYTES, MAX_TABLE_BYTES = 32 << 20, 64 << 20
MEMORY_BYTES, RESERVE = 4 << 30, 20 << 30
POLICY = 'largest current absolute interval width; tie order D,E,U; midpoint integer split'
CHECKPOINT_FORMAT = 'changed-nodes-v25; pending_frontier omitted only from checkpoints'
base, tail, cert = bands.base, bands.tail, bands.provider


def source_pins():
    for module, expected in ((bands, FROZEN_BANDS), (portable, FROZEN_PORTABLE)):
        if sha256(Path(module.__file__).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen v25 dependency changed: '+Path(module.__file__).name)
    pins = portable.source_pins()
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    for name, expected in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Pinned v25 source changed: '+name)


def check_pins(pins, witness_source):
    check_source_pins(pins)
    if sha256((HERE/witness_source['path']).read_bytes()).hexdigest() != witness_source['sha256']:
        raise RuntimeError('Pinned portable seed changed')


def split_axis(cell):
    """Use current widths, not ratios to the root; deterministic D/E/U ties."""
    if type(cell) is not bands.Rectangle:
        raise ValueError('Expected the frozen exact integer Rectangle')
    return max((1, 2, 0), key=lambda axis: cell.ranges[axis][1]-cell.ranges[axis][0])


def cover(evaluator, *, max_nodes=MAX_NODES, max_seconds=MAX_SECONDS,
          cell_failure_bits=72, aggregate_bits=54, checkpoint=None, context_checkpoint=None):
    """New controller, unchanged cell receipt and directed union operations.

    Injected evaluators are restricted to synthetic tests by the directed flag.
    A complete numerical band bound is emitted only after the frozen tree audit
    accepts the full integer partition. I/O/resource exceptions leave only
    partial checkpoint records, never a completed summary.
    """
    if (type(max_nodes) is not int or not 1 <= max_nodes <= MAX_NODES
            or type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= MAX_SECONDS
            or any(type(x) is not int or not 0 <= x <= 4096 for x in (cell_failure_bits, aggregate_bits))
            or any(callback is not None and not callable(callback) for callback in (checkpoint, context_checkpoint))):
        raise ValueError('Invalid v25 cover limits or callbacks')
    if getattr(evaluator, '_active', False):
        raise RuntimeError('A weight-band cover is already active')
    pins = source_pins()
    check_source_pins(pins)
    evaluator.check_pins()
    budget = base.Budget(max_seconds)
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
                changed_nodes=deepcopy(changed), evaluated_nodes=evaluations,
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
                if bands._passes(receipt, cell, cell_failure_bits, context_hash):
                    nodes[index] = dict(node, status='accepted', receipt=receipt)
                    publish([nodes[index]], 'accepted-leaf')
                elif receipt['complete_cell'] is not True:
                    nodes[index] = dict(node, status='unresolved', receipt=receipt)
                    blocked.append(index)
                    publish([nodes[index]], 'unresolved-leaf')
                else:
                    axis = split_axis(cell)
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
            accepted, pending = bands.audit_tree(nodes, roots, cell_failure_bits, context_hash)
            complete = not pending
            result = dict(schema='gamma-midpoint-weight-band-cover-v1', context=context,
                context_sha256=context_hash, status='complete' if complete else 'partial',
                requested_weight_range=context['requested_weight_range'],
                requested_weights=context['requested_weights'], nodes=nodes,
                root_cells=[cell.descriptor() for cell in roots], accepted_leaves=accepted,
                pending_frontier=pending, evaluated_nodes=evaluations, tree_nodes=len(nodes),
                checkpoints=sequence, complete_requested_band=complete,
                covers_all_supports_in_band=complete, covers_all_matching_profiles=complete,
                every_integer_weight_prepared=True,
                directed=complete and type(evaluator) is bands.Evaluator,
                numerical_backend='actual-all-pair-bank' if type(evaluator) is bands.Evaluator else 'synthetic-test-only',
                whole_code_certificate=False, no_weight_interpolation=True,
                stop_reason=None if complete else stop_reason or 'unresolved-leaves',
                zero=None, log2_upper_rational=None, log2_upper_enclosure=None,
                meets_requested_band_target=False, requested_failure_bits=aggregate_bits,
                cell_failure_bits=cell_failure_bits, max_nodes=max_nodes, max_seconds=max_seconds,
                controller=dict(path=NAME, sha256=pins[NAME]), controller_source_sha256=pins,
                partition_policy=POLICY, checkpoint_format=CHECKPOINT_FORMAT,
                construction_changed=False)
            if complete:
                # Identical positive directed union to frozen bands.cover.
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
        try:
            evaluator.check_pins()
        finally:
            check_source_pins(pins)
    result.update(batch_source_verified_at_entry=True, batch_source_verified_at_exit=True,
                  elapsed_seconds=monotonic()-budget.started)
    return result


def validate_band(first, last, max_nodes):
    if (type(first) is not int or type(last) is not int or not 368 <= first <= last <= 383
            or last-first+1 > 64 or type(max_nodes) is not int or not 1 <= max_nodes <= MAX_NODES):
        raise ValueError('Expected assigned weights in 368..383 and at most 4096 evaluated nodes')


def load_seed(name, expected_sha256):
    name = portable.relative_name(name, suffix='.json')
    portable.digest_string(expected_sha256)
    raw = (HERE/name).read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Portable seed byte pin mismatch')
    seed = json.loads(raw)
    parameters = portable.validate_seed(seed)
    return seed, parameters, dict(path=name, sha256=expected_sha256)


def reserve_guard(directory, planned_bytes=0):
    if type(planned_bytes) is not int or planned_bytes < 0:
        raise ValueError('Expected a nonnegative planned write size')
    if shutil.disk_usage(directory).free-planned_bytes < RESERVE:
        raise RuntimeError('Resource stop: planned output would leave fewer than 20 GiB free')


class OutputBudget:
    """Exact cumulative JSON bytes; exclusive writes, charged before opening."""

    def __init__(self, directory):
        self.directory, self.used_bytes = Path(directory).resolve(), 0

    def write_new(self, path, value):
        path = Path(path)
        if path.parent.resolve() != self.directory:
            raise ValueError('Artifact outside the new run directory')
        if path.exists():
            raise FileExistsError(path)
        size = len((json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8'))
        if size > MAX_OUTPUT_BYTES-self.used_bytes:
            raise RuntimeError('Resource stop: 32 MiB artifact budget; partial evidence preserved')
        reserve_guard(path.parent, size)
        self.used_bytes += size
        bands.profiles.write_new(path, value)


def install_memory_limit():
    """Fail closed off Linux; RLIMIT_AS also bounds native library allocations."""
    if sys.platform != 'linux':
        raise RuntimeError('Production v25 requires Linux RLIMIT_AS; pure tests are portable')
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    limit = min([MEMORY_BYTES]+[x for x in (soft, hard) if x != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (limit, hard))
    observed, unused = resource.getrlimit(resource.RLIMIT_AS)
    if observed == resource.RLIM_INFINITY or observed > MEMORY_BYTES:
        raise RuntimeError('Could not enforce the 4 GiB address-space limit')
    return dict(mechanism='Linux RLIMIT_AS', requested_bytes=MEMORY_BYTES, applied_bytes=observed)


def remaining_seconds(budget, ceiling=MAX_SECONDS):
    budget.check()
    return min(ceiling, budget.max_seconds-(monotonic()-budget.started))


def driver_record(pins, witness_source, first, last, max_nodes, seed_weight):
    return dict(schema='gamma-portable-compact-band-driver-v25',
        driver=dict(path=NAME, sha256=pins[NAME]), source_sha256=pins,
        witness_source=witness_source, source_calibration_weight=seed_weight,
        requested_weight_range=[first, last], max_nodes=max_nodes,
        max_total_compute_seconds=MAX_SECONDS, time_budget_kind='shared cooperative build/prepare/cover wall budget',
        max_table_bytes=MAX_TABLE_BYTES, max_address_space_bytes=MEMORY_BYTES,
        precision=256, coefficient_bits=80, grid=128, field_order=1 << 128,
        state_dimension=30, target_distance='2/5', construction=cert.local.construction(30),
        outer_measure='invertible-toeplitz', cell_failure_bits=72, aggregate_failure_bits=54,
        disk_reserve_bytes=RESERVE, max_output_bytes=MAX_OUTPUT_BYTES,
        output_budget_scope='Exact UTF-8 exclusive JSON writes; external console logs excluded',
        numerical_provider=Path(bands.__file__).name, partition_controller=NAME,
        partition_policy=POLICY, checkpoint_format=CHECKPOINT_FORMAT, final_frontier_preserved=True,
        fresh_bank_required=True, fresh_bank_rebuilt=False, old_numerical_bounds_reused=False,
        construction_changed=False, directed=False, whole_code_certificate=False)


def run(seed_name, expected_sha256, first, last, output_dir, *, max_nodes=MAX_NODES):
    validate_band(first, last, max_nodes)
    destination = HERE/portable.relative_name(output_dir)
    seed, (z, alpha, markers), witness_source = load_seed(seed_name, expected_sha256)
    pins = source_pins()
    check_pins(pins, witness_source)
    bank = None
    try:
        reserve_guard(HERE, MAX_OUTPUT_BYTES)
        memory = install_memory_limit()
        budget = base.Budget(MAX_SECONDS)
        destination.mkdir(exist_ok=False)
        artifacts = OutputBudget(destination)
        artifacts.write_new(destination/'driver.json', driver_record(
            pins, witness_source, first, last, max_nodes, seed['active_groups']))
        progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
        bank = cert.AliasBank(z, alpha, markers, state=30, q=1 << 128,
            distance=F(2, 5), grid=128, precision=256, bits=80,
            outer_measure='invertible-toeplitz', max_seconds=remaining_seconds(budget),
            max_table_bytes=MAX_TABLE_BYTES, progress=progress)
        check_pins(pins, witness_source)
        bank.check_pins()
        artifacts.write_new(destination/'bank.json', dict(bank.metadata(), witness_source=witness_source,
            runner_source_sha256=pins, fresh_bank_rebuilt=True, runtime=portable.runtime_metadata(),
            memory_limit=memory, directed=False, whole_code_certificate=False))
        evaluator = bands.Evaluator(bank, first, last, max_seconds=remaining_seconds(budget, 180.),
            witness_source=witness_source, progress=progress)
        result = cover(evaluator, max_nodes=max_nodes, max_seconds=remaining_seconds(budget),
            cell_failure_bits=72, aggregate_bits=54,
            checkpoint=lambda row: artifacts.write_new(destination/f"node_{row['sequence']:06d}.json", row),
            context_checkpoint=lambda row: artifacts.write_new(destination/'context.json', row))
        check_pins(pins, witness_source)
        bank.check_pins()
        result['total_compute_seconds'] = monotonic()-budget.started
        artifacts.write_new(destination/'summary.json', result)
        return result
    finally:
        try:
            check_pins(pins, witness_source)
        finally:
            if bank is not None:
                bank.check_pins()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--max-nodes', type=int, default=MAX_NODES)
    args = parser.parse_args()
    result = run(args.seed, args.expected_sha256, args.first, args.last,
                 args.output_dir, max_nodes=args.max_nodes)
    print(json.dumps({key: result[key] for key in ('status', 'requested_weight_range',
        'evaluated_nodes', 'log2_upper_rational', 'stop_reason')}, allow_nan=False), flush=True)
    return 0 if result['complete_requested_band'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
