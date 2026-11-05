"""Portable single-weight orchestration; frozen numerical producer and rows.

Only exact seed literals are used. Historical links are not followed.
Saved-data audit is separate; no whole-code claim.
"""
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
BASE = 'hd_field128_gamma30_d40_r913_midpoint_directed_v25'
BASE_SHA = '5c0f188939194cd5f9304644d25a11c8f0b1e97b57a0ce8eddc15f5157751a8b'
PORTABLE = 'hd_field128_gamma_portable_band_runner_v24'
PORTABLE_SHA = 'e69c97a06032bbe06aade7d9ed707f9353b633c464e70d740d8647ddea58ed1b'
MAX_SECONDS, MAX_OUTPUT_BYTES = 600, 2 << 20
MAX_TABLE_BYTES, MEMORY_BYTES, RESERVE = 64 << 20, 4 << 30, 20 << 30
Q, FAILURE_BITS, CUTOFF = 1 << 128, 54, 64
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'


def baseline():
    if sha256((HERE/(BASE+'.py')).read_bytes()).hexdigest() != BASE_SHA:
        raise RuntimeError('Frozen direct midpoint wrapper changed')
    return importlib.import_module(BASE)


def geometry(weight):
    if type(weight) is not int or not 1 <= weight <= 8191:
        raise ValueError('Expected one original weight in1..8191;8192 is a separate terminal event')
    first, last = max(0,weight-4096), weight//2
    rows = last-first+1
    return first, last, rows, rows*rows


def target_bits(value):
    if type(value) is not int or not 40 <= value <= 4096:
        raise ValueError('Expected an integer failure target in40..4096')
    return value


def source_pins(previous=None, portable=None, producer=None):
    previous = baseline() if previous is None else previous
    portable = previous.pinned_module(PORTABLE, PORTABLE_SHA) if portable is None else portable
    producer = previous.pinned_module(previous.PRODUCER, previous.PRODUCER_SHA256) if producer is None else producer
    # No source map or historical link from the seed enters this closure.
    pins = previous.source_pins(producer=producer)
    for name, digest in portable.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def driver_record(pins, witness_source, seed, weight, memory_limit, failure_bits=FAILURE_BITS):
    first, last, rows, profiles = geometry(weight)
    return dict(schema='gamma-portable-selected-midpoint-driver-v29',
        driver=dict(path=NAME,sha256=pins[NAME]), source_sha256=pins,
        witness_source=witness_source, active_groups=weight, requested_weight_range=[weight,weight],
        source_calibration_weight=seed['active_groups'], construction=CONSTRUCTION,
        field_order=Q, state_dimension=30, target_distance='2/5',
        marker_count=len(seed['witness']['markers']), numerical_provider='hd_field128_gamma_midpoint_certificate.py',
        expected_even_full_range=[first,last], expected_even_rows=rows, expected_profile_pairs=profiles,
        precision=256, coefficient_bits=80, grid=128, requested_failure_bits=target_bits(failure_bits),
        maximum_refined_singletons=CUTOFF, whole_job_seconds=MAX_SECONDS,
        max_table_bytes=MAX_TABLE_BYTES, max_output_bytes=MAX_OUTPUT_BYTES,
        disk_reserve_bytes=RESERVE, memory_address_space_bytes=memory_limit,
        required_normal_process_exit=True, fresh_bank_required=True, fresh_bank_rebuilt=False,
        literal_seed_only=True, historical_origin_verified=False, calibration_history_replayed=False,
        old_numerical_bounds_reused=False, directed=False, whole_code_certificate=False)


def admitted(result, weight, failure_bits=FAILURE_BITS):
    first, last, rows, profiles = geometry(weight)
    failure_bits = target_bits(failure_bits)
    required = dict(status='complete', failure=None, active_groups=weight, directed=True,
        complete_selected_weight=True, covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
        whole_code_certificate=False, completed_even_rows=rows, expected_even_rows=rows,
        completed_even_full_range=[first,last], pending_even_full_range=None,
        profile_pairs_evaluated=profiles, expected_profile_pairs=profiles,
        requested_failure_bits=failure_bits, meets_requested_selected_weight_target=True)
    if any(key not in result or result[key] != value for key,value in required.items()):
        return False
    if any(type(result.get(key)) is not type(value) for key,value in required.items()):
        return False
    try:
        return F(result['log2_upper_rational']) <= -failure_bits
    except (KeyError,TypeError,ValueError,ZeroDivisionError):
        return False


def run(seed_name, expected_sha256, weight, output_dir, *, failure_bits=FAILURE_BITS):
    geometry(weight)
    failure_bits = target_bits(failure_bits)
    previous = baseline()  # Standard-library-only module; no heavy imports yet.
    limits = previous.validation()
    seed_name = limits.relative_name(seed_name,json_file=True,basename=True)
    limits.digest(expected_sha256)
    destination = HERE/limits.relative_name(output_dir,basename=True)
    limits.reserve_guard(MAX_OUTPUT_BYTES)
    started, bank = monotonic(), None
    with limits.linux_limits() as memory_limit:
        portable = previous.pinned_module(PORTABLE,PORTABLE_SHA)
        producer = previous.pinned_module(previous.PRODUCER,previous.PRODUCER_SHA256)
        seed, (z,alpha,markers), witness_source = portable.load_seed(seed_name,expected_sha256)
        pins = source_pins(previous,portable,producer)
        limits.check_pins(pins,witness_source)
        destination.mkdir(exist_ok=False)
        outputs = previous.OutputBudget(destination,limits.reserve_guard,producer.write_new)
        written = []

        def save(name,value,*,phase_check=True):
            if phase_check:
                limits.check_pins(pins,witness_source)
            path = destination/name
            outputs.write_new(path,value)
            written.append((path,sha256(path.read_bytes()).hexdigest()))

        def remaining():
            seconds = MAX_SECONDS-5-(monotonic()-started)
            if seconds <= 0:
                raise TimeoutError('Whole-job budget exhausted')
            return seconds

        try:
            save('driver.json',driver_record(pins,witness_source,seed,weight,memory_limit,failure_bits))
            bank = producer.provider.AliasBank(z,alpha,markers,state=30,q=Q,distance=F(2,5),
                grid=128,precision=256,bits=80,outer_measure='invertible-toeplitz',
                max_seconds=remaining(),max_table_bytes=MAX_TABLE_BYTES)
            bank.check_pins()
            save('bank.json',dict(bank.metadata(),witness_source=witness_source,runner_source_sha256=pins,
                runtime=portable.runtime_metadata(),fresh_bank_rebuilt=True,directed=False,whole_code_certificate=False))
            result = producer.weight_receipt(bank,weight,max_singletons=CUTOFF,
                max_seconds=remaining(),failure_bits=failure_bits,witness_source=witness_source,
                checkpoint=lambda row: save(f"row_{row['even_full']:04d}.json",row,phase_check=False),
                context_checkpoint=lambda row: save('context.json',row))
            bank.check_pins()
            result.update(execution_driver=dict(path=NAME,sha256=pins[NAME]),execution_source_sha256=pins,
                source_calibration_weight=seed['active_groups'],literal_seed_only=True,
                historical_origin_verified=False,calibration_history_replayed=False,
                required_normal_process_exit=True,wrapper_admitted=admitted(result,weight,failure_bits),
                total_job_seconds=monotonic()-started)
            save('summary.json',result)
            return result
        finally:
            try:
                limits.check_pins(pins,witness_source)
                for path,expected in written:
                    if sha256(path.read_bytes()).hexdigest() != expected:
                        raise RuntimeError('New output bytes changed: '+path.name)
            finally:
                if bank is not None:
                    bank.check_pins()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',required=True)
    parser.add_argument('--expected-sha256',required=True)
    parser.add_argument('--weight',type=int,required=True)
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--failure-bits',type=int,default=FAILURE_BITS)
    args = parser.parse_args()
    result = run(args.seed,args.expected_sha256,args.weight,args.output_dir,failure_bits=args.failure_bits)
    print(json.dumps({key:result.get(key) for key in ('status','active_groups','completed_even_rows',
        'profile_pairs_evaluated','log2_upper_rational','meets_requested_selected_weight_target','wrapper_admitted')}),flush=True)
    return 0 if admitted(result,args.weight,args.failure_bits) else 2


if __name__ == '__main__':
    raise SystemExit(main())
