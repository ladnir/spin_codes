"""Portable bands with an explicit cell threshold and a fixed 54-bit union gate.

The frozen v26 scope checks and v25 controller/evaluator are unchanged. The
cell threshold controls acceptance and splitting, not the numerical bound.
A complete tree is successful only when its saved whole-band union meets
54 bits. No inherited numerical bound or whole-code claim is imported.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from time import monotonic

import hd_field128_gamma_portable_band_runner_v26 as previous


controller = previous.controller
HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FROZEN_PREVIOUS = 'c68765dad2b5fbca7e8fe0cd3a1d11c977695c243ac3951a72f5a405f296b26e'
AGGREGATE_BITS = 54


def validate(first, last, max_nodes, cell_failure_bits):
    previous.validate_band(first, last, max_nodes)
    if type(cell_failure_bits) is not int or not 54 <= cell_failure_bits <= 128:
        raise ValueError('Expected integer cell_failure_bits in 54..128')


def source_pins():
    if sha256(Path(previous.__file__).read_bytes()).hexdigest() != FROZEN_PREVIOUS:
        raise RuntimeError('Frozen v26 portable runner changed')
    pins = previous.source_pins()
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


@contextmanager
def hard_deadline():
    """Hard process exit cannot be swallowed by a cooperative timeout handler."""
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Actual portable runs require Linux; pure tests mock this context')
    old_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, lambda unused_signum, unused_frame: os._exit(124))
    signal.setitimer(signal.ITIMER_REAL, controller.MAX_SECONDS)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)


def driver_record(pins, witness_source, first, last, max_nodes, seed_weight, cell_failure_bits=72):
    validate(first, last, max_nodes, cell_failure_bits)
    record = previous.driver_record(pins, witness_source, first, last, max_nodes, seed_weight)
    record.update(schema='gamma-portable-compact-band-driver-v27',
        driver=dict(path=NAME, sha256=pins[NAME]), source_sha256=pins,
        cell_failure_bits=cell_failure_bits, aggregate_failure_bits=AGGREGATE_BITS,
        allowed_cell_failure_bits=[54, 128], hard_wall_seconds=controller.MAX_SECONDS,
        success_requires_complete_band_and_54_bit_union=True,
        wrapper_scope='Cell acceptance threshold only; frozen v25 cover and numerical bounds unchanged')
    return record


def admitted(result):
    """Check actual aggregate admission, not just a complete partition flag."""
    if (result.get('status') != 'complete' or result.get('complete_requested_band') is not True
            or result.get('directed') is not True or result.get('whole_code_certificate') is not False
            or result.get('requested_failure_bits') != AGGREGATE_BITS
            or result.get('meets_requested_band_target') is not True
            or result.get('pending_frontier') != []):
        return False
    if result.get('zero') is True:
        return result.get('log2_upper_rational') is None
    if result.get('zero', False) is not False:
        return False
    try:
        return F(result['log2_upper_rational']) <= -AGGREGATE_BITS
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def run(seed_name, expected_sha256, first, last, output_dir, *,
        max_nodes=controller.MAX_NODES, cell_failure_bits=72):
    validate(first, last, max_nodes, cell_failure_bits)
    with hard_deadline():
        destination = HERE/controller.portable.relative_name(output_dir)
        seed, (z, alpha, markers), witness_source = controller.load_seed(seed_name, expected_sha256)
        pins = source_pins()
        previous.check_pins(pins, witness_source)
        bank = None
        try:
            controller.reserve_guard(HERE, controller.MAX_OUTPUT_BYTES)
            memory = controller.install_memory_limit()
            budget = controller.base.Budget(controller.MAX_SECONDS)
            destination.mkdir(exist_ok=False)
            artifacts = controller.OutputBudget(destination)
            artifacts.write_new(destination/'driver.json', driver_record(
                pins, witness_source, first, last, max_nodes, seed['active_groups'], cell_failure_bits))
            progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
            bank = controller.cert.AliasBank(z, alpha, markers, state=30, q=1 << 128,
                distance=F(2, 5), grid=128, precision=256, bits=80,
                outer_measure='invertible-toeplitz', max_seconds=controller.remaining_seconds(budget),
                max_table_bytes=controller.MAX_TABLE_BYTES, progress=progress)
            previous.check_pins(pins, witness_source)
            bank.check_pins()
            artifacts.write_new(destination/'bank.json', dict(bank.metadata(), witness_source=witness_source,
                runner_source_sha256=pins, fresh_bank_rebuilt=True,
                runtime=controller.portable.runtime_metadata(), memory_limit=memory,
                directed=False, whole_code_certificate=False))
            evaluator = controller.bands.Evaluator(bank, first, last,
                max_seconds=controller.remaining_seconds(budget, 180.),
                witness_source=witness_source, progress=progress)
            result = controller.cover(evaluator, max_nodes=max_nodes,
                max_seconds=controller.remaining_seconds(budget),
                cell_failure_bits=cell_failure_bits, aggregate_bits=AGGREGATE_BITS,
                checkpoint=lambda row: artifacts.write_new(destination/f"node_{row['sequence']:06d}.json", row),
                context_checkpoint=lambda row: artifacts.write_new(destination/'context.json', row))
            previous.check_pins(pins, witness_source)
            bank.check_pins()
            result.update(runner=dict(path=NAME, sha256=pins[NAME]), runner_source_sha256=pins,
                runner_admitted=admitted(result), total_compute_seconds=monotonic()-budget.started)
            artifacts.write_new(destination/'summary.json', result)
            return result
        finally:
            try:
                previous.check_pins(pins, witness_source)
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
    parser.add_argument('--cell-failure-bits', type=int, default=72)
    args = parser.parse_args()
    result = run(args.seed, args.expected_sha256, args.first, args.last, args.output_dir,
        max_nodes=args.max_nodes, cell_failure_bits=args.cell_failure_bits)
    print(json.dumps({key: result[key] for key in ('status', 'requested_weight_range',
        'evaluated_nodes', 'log2_upper_rational', 'stop_reason', 'meets_requested_band_target')},
        allow_nan=False), flush=True)
    return 0 if admitted(result) else 2


if __name__ == '__main__':
    raise SystemExit(main())
