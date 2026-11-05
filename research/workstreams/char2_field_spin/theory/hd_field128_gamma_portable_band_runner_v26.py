"""Portable explicit-band orchestration using the frozen v25 controller.

Accept 1..64 consecutive original weights in 1..8191. Weight 8192 is not
supported by the frozen regional evaluator and belongs to a separate terminal
event; this runner never silently omits it. No construction, coefficient,
cell-bound, splitting, checkpoint or union arithmetic changes in this version.
Only exact seed literals are imported, never an ancestral numerical bound.

Unchanged limits: 4096 evaluated nodes, shared cooperative 600-second compute
budget, Linux 4 GiB address space, 64 MiB coefficient tables, 32 MiB JSON output,
and 20 GiB remaining disk after each planned write. Root coordinates execution.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma_portable_band_runner_v25 as controller


HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FROZEN_CONTROLLER = '5d22c83e3419b1dc2a484b09c9cb2ada17971bd48195987140c9158942d7ef1c'


def validate_band(first, last, max_nodes):
    if (type(first) is not int or type(last) is not int
            or not 1 <= first <= last <= 8191 or last-first+1 > 64
            or type(max_nodes) is not int or not 1 <= max_nodes <= controller.MAX_NODES):
        raise ValueError('Expected 1..64 consecutive weights in 1..8191 and 1..4096 nodes; '
                         '8192 requires the separate terminal event')


def source_pins():
    if sha256(Path(controller.__file__).read_bytes()).hexdigest() != FROZEN_CONTROLLER:
        raise RuntimeError('Frozen v25 controller changed')
    pins = controller.source_pins()
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins, witness_source):
    # Both files live next to the frozen controller; all emitted paths remain
    # relative to that theory directory on the remote host and after transfer.
    controller.check_pins(pins, witness_source)


def driver_record(pins, witness_source, first, last, max_nodes, seed_weight):
    validate_band(first, last, max_nodes)
    record = controller.driver_record(pins, witness_source, first, last, max_nodes, seed_weight)
    record.update(schema='gamma-portable-compact-band-driver-v26',
        driver=dict(path=NAME, sha256=pins[NAME]), source_sha256=pins,
        partition_controller=controller.NAME,
        partition_controller_sha256=pins[controller.NAME],
        allowed_weight_range=[1, 8191], maximum_band_width=64,
        terminal_weight_8192_supported=False,
        wrapper_scope='Explicit band orchestration only; frozen v25 controller and numerical provider unchanged')
    return record


def run(seed_name, expected_sha256, first, last, output_dir, *, max_nodes=controller.MAX_NODES):
    validate_band(first, last, max_nodes)
    destination = HERE/controller.portable.relative_name(output_dir)
    seed, (z, alpha, markers), witness_source = controller.load_seed(seed_name, expected_sha256)
    pins = source_pins()
    check_pins(pins, witness_source)
    bank = None
    try:
        controller.reserve_guard(HERE, controller.MAX_OUTPUT_BYTES)
        memory = controller.install_memory_limit()
        budget = controller.base.Budget(controller.MAX_SECONDS)
        destination.mkdir(exist_ok=False)
        artifacts = controller.OutputBudget(destination)
        artifacts.write_new(destination/'driver.json', driver_record(
            pins, witness_source, first, last, max_nodes, seed['active_groups']))
        progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
        bank = controller.cert.AliasBank(z, alpha, markers, state=30, q=1 << 128,
            distance=F(2, 5), grid=128, precision=256, bits=80,
            outer_measure='invertible-toeplitz', max_seconds=controller.remaining_seconds(budget),
            max_table_bytes=controller.MAX_TABLE_BYTES, progress=progress)
        check_pins(pins, witness_source)
        bank.check_pins()
        artifacts.write_new(destination/'bank.json', dict(bank.metadata(), witness_source=witness_source,
            runner_source_sha256=pins, fresh_bank_rebuilt=True,
            runtime=controller.portable.runtime_metadata(), memory_limit=memory,
            directed=False, whole_code_certificate=False))
        evaluator = controller.bands.Evaluator(bank, first, last,
            max_seconds=controller.remaining_seconds(budget, 180.),
            witness_source=witness_source, progress=progress)
        result = controller.cover(evaluator, max_nodes=max_nodes,
            max_seconds=controller.remaining_seconds(budget), cell_failure_bits=72, aggregate_bits=54,
            checkpoint=lambda row: artifacts.write_new(destination/f"node_{row['sequence']:06d}.json", row),
            context_checkpoint=lambda row: artifacts.write_new(destination/'context.json', row))
        check_pins(pins, witness_source)
        bank.check_pins()
        result.update(runner=dict(path=NAME, sha256=pins[NAME]), runner_source_sha256=pins,
                      total_compute_seconds=monotonic()-budget.started)
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
    parser.add_argument('--max-nodes', type=int, default=controller.MAX_NODES)
    args = parser.parse_args()
    result = run(args.seed, args.expected_sha256, args.first, args.last,
                 args.output_dir, max_nodes=args.max_nodes)
    print(json.dumps({key: result[key] for key in ('status', 'requested_weight_range',
        'evaluated_nodes', 'log2_upper_rational', 'stop_reason')}, allow_nan=False), flush=True)
    return 0 if result['complete_requested_band'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
