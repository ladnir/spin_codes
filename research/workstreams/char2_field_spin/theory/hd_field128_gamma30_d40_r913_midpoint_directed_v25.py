"""Bounded original-weight913 all-profile replay with frozen midpoint counts."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib
import json
from pathlib import Path
from time import monotonic

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
VALIDATOR = 'hd_field128_gamma30_d40_r913_directed_v25'
VALIDATOR_SHA256 = 'e9302f55a676f6ed4a0efd1a032a103e6a020f1099456ab40a79d70ae0731e1d'
PRODUCER = 'hd_field128_gamma_midpoint_certificate'
PRODUCER_SHA256 = 'ff8b6bb70597dc06ec678d7b9b1dda935e0e2c6e64d15f25d4f1767fb9a0f8cb'
MAX_OUTPUT_BYTES = 2 << 20
EXPECTED_ROWS, EXPECTED_PROFILES = 457, 208849


def pinned_module(name, expected):
    if sha256((HERE/(name+'.py')).read_bytes()).hexdigest() != expected:
        raise RuntimeError('Reviewed helper changed: '+name)
    return importlib.import_module(name)


def validation():
    # This module imports only the standard library; heavy helpers load inside
    # the inherited hard600s/4GiB context, not before memory limits are installed.
    return pinned_module(VALIDATOR, VALIDATOR_SHA256)


def source_pins(proposal=None, validator=None, producer=None):
    validator = validation() if validator is None else validator
    producer = pinned_module(PRODUCER, PRODUCER_SHA256) if producer is None else producer
    pins = validator.source_pins(proposal)
    for name, digest in producer.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


class OutputBudget:
    """Exact aggregate bytes, including duplicated final row records."""
    def __init__(self, directory, reserve_guard, writer):
        self.directory = Path(directory).resolve()
        self.reserve_guard, self.writer, self.used_bytes = reserve_guard, writer, 0

    def write_new(self, path, value):
        path = Path(path)
        if path.parent.resolve() != self.directory:
            raise ValueError('Artifact outside new run directory')
        if path.exists():
            raise FileExistsError(path)
        size = len((json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8'))
        if size > MAX_OUTPUT_BYTES-self.used_bytes:
            raise RuntimeError('Resource stop:2MiB output budget; no complete summary')
        self.reserve_guard(size)
        self.used_bytes += size
        self.writer(path, value)


def admitted(result):
    """Never accept a partial receipt or merely a producer target flag."""
    return (result.get('status') == 'complete'
        and result.get('complete_selected_weight') is True
        and result.get('covers_all_supports_at_weight') is True
        and result.get('covers_all_even_odd_splits') is True
        and result.get('directed') is True
        and result.get('whole_code_certificate') is False
        and result.get('active_groups') == 913
        and result.get('completed_even_rows') == EXPECTED_ROWS
        and result.get('profile_pairs_evaluated') == EXPECTED_PROFILES
        and result.get('pending_even_full_range') is None
        and result.get('requested_failure_bits') == 54
        and result.get('meets_requested_selected_weight_target') is True
        and result.get('log2_upper_rational') is not None
        and F(result['log2_upper_rational']) <= -54)


def run(proposal_name, expected_sha256, output_dir):
    validator = validation()
    destination = HERE/validator.relative_name(output_dir, basename=True)
    proposal, (z, alpha, markers), witness_source = validator.load_proposal(proposal_name, expected_sha256)
    validator.reserve_guard(MAX_OUTPUT_BYTES)
    started, bank = monotonic(), None
    with validator.linux_limits() as memory_limit:
        producer = pinned_module(PRODUCER, PRODUCER_SHA256)
        pins = source_pins(proposal, validator, producer)
        validator.check_pins(pins, witness_source)
        destination.mkdir(exist_ok=False)
        outputs = OutputBudget(destination, validator.reserve_guard, producer.write_new)
        written = []

        def save(name, value, *, phase_check=True):
            if phase_check:
                validator.check_pins(pins, witness_source)
            path = destination/name
            outputs.write_new(path, value)
            written.append((path, sha256(path.read_bytes()).hexdigest()))

        def remaining():
            seconds = validator.MAX_SECONDS-5-(monotonic()-started)
            if seconds <= 0:
                raise TimeoutError('Whole-job budget exhausted')
            return seconds

        try:
            save('driver.json', dict(schema='gamma-r913-midpoint-directed-driver-v25',
                driver=dict(path=NAME, sha256=pins[NAME]), source_sha256=pins,
                witness_source=witness_source, active_groups=913, construction=validator.CONSTRUCTION,
                field_order=validator.Q, state_dimension=30, target_distance='2/5', marker_count=83,
                numerical_provider=PRODUCER+'.py', expected_even_rows=EXPECTED_ROWS,
                expected_profile_pairs=EXPECTED_PROFILES, precision=256, coefficient_bits=80,
                grid=128, requested_failure_bits=54, maximum_refined_singletons=64,
                whole_job_seconds=validator.MAX_SECONDS, max_table_bytes=validator.MAX_TABLE_BYTES,
                max_output_bytes=MAX_OUTPUT_BYTES, disk_reserve_bytes=validator.RESERVE,
                memory_address_space_bytes=memory_limit, required_normal_process_exit=True,
                fresh_bank_required=True, fresh_bank_rebuilt=False, old_numerical_bounds_reused=False,
                directed=False, whole_code_certificate=False))
            bank = producer.provider.AliasBank(z, alpha, markers, state=30, q=validator.Q,
                distance=F(2,5), grid=128, precision=256, bits=80,
                outer_measure='invertible-toeplitz', max_seconds=remaining(),
                max_table_bytes=validator.MAX_TABLE_BYTES)
            bank.check_pins()
            save('bank.json', dict(bank.metadata(), witness_source=witness_source,
                runner_source_sha256=pins, runtime=validator.helper().portable.runtime_metadata(),
                fresh_bank_rebuilt=True, directed=False, whole_code_certificate=False))
            result = producer.weight_receipt(bank, 913, max_singletons=64,
                max_seconds=remaining(), failure_bits=54, witness_source=witness_source,
                # The frozen producer already checks its numerical closure at
                # every row; wrapper closure checks bracket this entire call.
                checkpoint=lambda row: save(f"row_{row['even_full']:04d}.json", row, phase_check=False),
                context_checkpoint=lambda row: save('context.json', row))
            bank.check_pins()
            result.update(execution_driver=dict(path=NAME, sha256=pins[NAME]),
                execution_source_sha256=pins, required_normal_process_exit=True,
                total_job_seconds=monotonic()-started)
            save('summary.json', result)
            return result
        finally:
            try:
                validator.check_pins(pins, witness_source)
                for path, expected in written:
                    if sha256(path.read_bytes()).hexdigest() != expected:
                        raise RuntimeError('New output bytes changed: '+path.name)
            finally:
                if bank is not None:
                    bank.check_pins()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    result = run(args.proposal, args.expected_sha256, args.output_dir)
    print(json.dumps({key: result.get(key) for key in ('status', 'active_groups',
        'completed_even_rows', 'profile_pairs_evaluated', 'log2_upper_rational',
        'meets_requested_selected_weight_target')}), flush=True)
    return 0 if admitted(result) else 2


if __name__ == '__main__':
    raise SystemExit(main())
