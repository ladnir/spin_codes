"""Audit fresh portable band receipts, not their seed's calibration history.

This distinct audit identity trusts only the byte-pinned positive seed literals
and reviewed fresh numerical producer. Frozen metadata, count, partition and
saved-enclosure union checks are reused. No coefficients or transfers replay.
Historical seed ancestry is explicitly outside scope and is never opened.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_midpoint_weight_band_typed_v2_audit as typed
import hd_field128_gamma_portable_band_runner_v27 as runner

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
TYPED_SHA256 = '5275bfd968cfdfdf79952137fdd787e1cdd205bc6cb5b977887142e7774a9095'
RUNNER_SHA256 = '87746a1d6a54f01baf844a3f9b3d02c4963fad0a26895f0051ab5711a0e26a1f'
SCHEMA = 'gamma-midpoint-weight-band-literal-only-remote-audit-v1'
MAX_OUTPUT_BYTES, RESERVE = 1 << 20, 20 << 30
require = typed.require


def source_pins():
    for module, digest in ((typed, TYPED_SHA256), (runner, RUNNER_SHA256)):
        require(sha256(Path(module.__file__).read_bytes()).hexdigest() == digest,
                'Frozen audit/runner identity changed')
    pins = typed.source_pins()
    for name, digest in runner.source_pins().items():
        require(name not in pins or pins[name] == digest, 'Conflicting audit source pin')
        pins[name] = digest
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def _sources(value, pins):
    """Validate fresh source maps only; never interpret historical path links."""
    if isinstance(value, dict):
        for name, child in value.items():
            if name == 'source_sha256' or name.endswith('_source_sha256'):
                typed._source_map(child, pins)
            else:
                _sources(child, pins)
    elif isinstance(value, list):
        for child in value:
            _sources(child, pins)


def _checkpoints(record, directory, inputs, budget):
    # Exact frozen typed-v2 checkpoint reconstruction, including partial flags.
    count = record['checkpoints']
    require(type(count) is int and count == len(record['nodes'])+1, 'Missing node checkpoints')
    replay = {}
    for sequence in range(count):
        budget.check()
        checkpoint, unused = inputs.load(directory/f'node_{sequence:06d}.json')
        require(checkpoint['sequence'] == sequence
                and checkpoint['context_sha256'] == record['context_sha256'], 'Checkpoint linkage mismatch')
        for key in ('complete_requested_band', 'covers_all_supports_in_band', 'directed',
                    'batch_source_verified_at_exit', 'whole_code_certificate'):
            require(checkpoint[key] is False, 'Checkpoint overclaims completion')
        require(checkpoint['log2_upper_rational'] is None, 'Partial checkpoint has band bound')
        for row in checkpoint['changed_nodes']:
            replay[row['id']] = row
    require([replay[i] for i in sorted(replay)] == record['nodes'], 'Checkpoint replay differs from final tree')
    return count


def _seed(context, inputs, expected_sha256):
    portable = runner.controller.portable
    portable.digest_string(expected_sha256)
    pin = context.get('witness_source')
    require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'Missing exact portable seed pin')
    require(pin['sha256'] == expected_sha256, 'Portable seed differs from explicit expected digest')
    name = portable.relative_name(pin['path'], suffix='.json')
    seed, seed_pin = inputs.load(HERE/name, expected_sha256)
    portable.validate_seed(seed)
    require(seed.get('selected_weight_certificate') is False, 'Seed cannot claim a certificate')
    for key in ('z', 'alpha', 'markers', 'grid', 'state_dimension'):
        require(seed['witness'][key] == context['witness'][key], 'Portable seed/context literal mismatch: '+key)
    # No _provenance(seed), source-map recursion, or historical-origin reads.
    return seed, seed_pin


def audit_summary(path, *, expected_summary_sha256, expected_seed_sha256, max_seconds=120.):
    portable = runner.controller.portable
    portable.digest_string(expected_summary_sha256)
    budget, pins, inputs = typed.base.Budget(max_seconds), source_pins(), typed.Inputs()
    typed.check_source_pins(pins)
    record, summary_pin = inputs.load(path, expected_summary_sha256)
    typed.previous._metadata(record)
    require(runner.admitted(record), 'Only complete 54-bit admitted bands can be audited')
    context = record['context']
    seed, seed_pin = _seed(context, inputs, expected_seed_sha256)
    directory = Path(summary_pin['path']).parent
    driver, driver_pin = inputs.load(directory/'driver.json')
    bank, bank_pin = inputs.load(directory/'bank.json')
    runner_pins = runner.source_pins()
    first, last = record['requested_weight_range']
    expected_driver = runner.driver_record(runner_pins, context['witness_source'], first, last,
        record['max_nodes'], seed['active_groups'], record['cell_failure_bits'])
    require(driver == expected_driver, 'Exact v27 execution driver differs')
    require(record['runner'] == driver['driver'] and record['runner_source_sha256'] == runner_pins,
            'Summary runner identity differs')
    require(record['controller'] == dict(path=runner.controller.NAME,
            sha256=runner_pins[runner.controller.NAME]), 'Controller identity differs')
    require(record['controller_source_sha256'] == runner.controller.source_pins(),
            'Controller source closure differs')
    require(record['runner_admitted'] is True and record['construction_changed'] is False,
            'Runner admission/construction mismatch')
    require(record['total_compute_seconds'] <= driver['max_total_compute_seconds'], 'Reported compute exceeds limit')
    require(bank['runner_source_sha256'] == runner_pins
            and bank['witness_source'] == context['witness_source'], 'Bank runner/witness pin differs')
    require(bank['fresh_bank_rebuilt'] is True and bank['directed'] is False
            and bank['whole_code_certificate'] is False, 'Wrong fresh bank metadata scope')
    require(bank['witness'] == context['witness'], 'Bank/context exact witness differs')
    for key in ('construction', 'field_order', 'state_dimension', 'target_distance', 'output_length',
                'source_dimension', 'outer_measure', 'bank_source_sha256', 'all_pair_before_alpha',
                'old_numerical_ceilings_reused', 'construction_unchanged'):
        require(bank[key] == context[key], 'Bank/context construction differs: '+key)
    require(bank['memory_limit'] == dict(mechanism='Linux RLIMIT_AS', requested_bytes=4<<30,
            applied_bytes=4<<30), 'Wrong reported memory limit')
    require(bank['runtime']['OPENBLAS_NUM_THREADS'] == '1', 'Wrong reported BLAS thread count')
    require(context['witness']['precision'] == 256 and context['witness']['grid'] == 128
            and context['witness']['packet_mantissa_bits'] == 80, 'Wrong fresh bank precision/grid')
    saved_context, context_pin = inputs.load(directory/'context.json')
    require(saved_context == dict(context=context, context_sha256=record['context_sha256'], directed=False,
            complete_requested_band=False, whole_code_certificate=False), 'Context artifact mismatch')
    count = _checkpoints(record, directory, inputs, budget)
    sources = dict(pins)
    for value in (context, driver, bank, {k:v for k,v in record.items() if k not in ('context','nodes')}):
        _sources(value, sources)
    # Exact original-r normalization, complete roots/tree, all stored counts,
    # accepted leaf enclosures and fresh 384-bit union: unchanged frozen code.
    numeric = typed.previous.audit_numeric(record, 4096, budget)
    artifact_paths = [directory/name for name in ('driver.json','bank.json','context.json','summary.json')]
    artifact_paths += [directory/f'node_{i:06d}.json' for i in range(count)]
    output_bytes = sum(p.stat().st_size for p in artifact_paths)
    require(output_bytes <= driver['max_output_bytes'], 'Producer JSON exceeded declared output budget')
    typed.check_source_pins(sources)
    inputs.check()
    return dict(schema=SCHEMA, audit_result='PASS', input_summary=summary_pin, context_file=context_pin,
        execution_driver_receipt=driver_pin, bank_metadata=bank_pin, portable_seed=seed_pin,
        context_sha256=record['context_sha256'], construction=typed.CONSTRUCTION,
        state_dimension=30, field_order=1<<128, target_distance='2/5',
        requested_weight_range=record['requested_weight_range'], requested_weights=numeric['weights'],
        complete_requested_band=True, covers_all_matching_profiles=True, integer_partition_checked=True,
        every_original_weight_preparation_metadata_checked=True, normalized_histograms=numeric['normalized_histograms'],
        counts_recomputed=True, count_algorithm='unchanged frozen independent midpoint prefix recurrence',
        cells_checked=numeric['cells_checked'], accepted_cells=numeric['accepted_cells'], checkpoints_checked=count,
        cell_log_enclosures_checked=True, positive_count_zero_cells=numeric['positive_count_zero_cells'],
        directed_union_precision=384, saved_band_log2_upper_rational=record['log2_upper_rational'],
        reaggregated_cell_log2_upper_rational=numeric['union'], reaggregated_cell_log2_upper_enclosure=numeric['union_enclosure'],
        log2_upper_rational=record['log2_upper_rational'], zero=numeric['zero'],
        source_sha256=pins, checked_source_sha256=dict(sorted(sources.items())),
        checked_file_pins=list(inputs.pins.values()), source_and_input_bytes_checked_at_exit=True,
        frozen_typed_v2_auditor_sha256=TYPED_SHA256, frozen_numerical_auditor_sha256=typed.PREVIOUS_SHA256,
        exact_seed_literals_checked=True, exact_driver_record_checked=True, fresh_bank_metadata_checked=True,
        producer_json_bytes=output_bytes, producer_json_budget_bytes=driver['max_output_bytes'],
        source_calibration_weight=seed['active_groups'], historical_origin_verified=False,
        calibration_history_replayed=False, historical_bounds_used=False,
        provenance_scope='Pinned portable seed is a literal-only input; historical ancestry is not a mathematical premise and is not verified',
        source_pin_conflicts=False, coefficient_replay_in_this_audit=False,
        normalized_regional_replay=False, cell_transfer_recomputed=False,
        trusted_numerical_producer='Reviewed byte-pinned fresh band producer; saved caps and resource telemetry are not authenticated by this audit',
        no_weight_interpolation=True, whole_code_certificate=False)


def write_remote(path, report):
    raw = (json.dumps(report, indent=2, allow_nan=False)+'\n').encode()
    require(len(raw) <= MAX_OUTPUT_BYTES, 'Audit output exceeds one MiB')
    require(shutil.disk_usage(Path(path).parent).free-len(raw) >= RESERVE, 'Audit output breaches 20 GiB reserve')
    return typed.producer.profiles.write_new(path, report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('--expected-summary-sha256', required=True)
    parser.add_argument('--expected-seed-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120.)
    args = parser.parse_args()
    require(not args.output.exists(), 'Audit output already exists')
    require(shutil.disk_usage(args.output.parent).free-MAX_OUTPUT_BYTES >= RESERVE,
            'Audit output preflight breaches 20 GiB reserve')
    result = audit_summary(args.summary, expected_summary_sha256=args.expected_summary_sha256,
        expected_seed_sha256=args.expected_seed_sha256, max_seconds=args.max_seconds)
    print(json.dumps(dict(artifact=write_remote(args.output, result),
        cells_checked=result['cells_checked'], accepted_cells=result['accepted_cells'],
        checkpoints_checked=result['checkpoints_checked'], log2_upper_rational=result['log2_upper_rational'],
        historical_origin_verified=False)), flush=True)
