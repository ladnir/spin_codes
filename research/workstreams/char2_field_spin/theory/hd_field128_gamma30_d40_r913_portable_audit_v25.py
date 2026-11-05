"""Audit the saved portable913 receipt without rewriting its relative links.

The frozen all-profile count and row/prefix checker is unchanged. This adapter
resolves portable proposal links from theory, as the producing wrapper does,
and independently binds the external context and all457 row checkpoints.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath, PureWindowsPath
from time import monotonic

import hd_field128_gamma_midpoint_receipt_audit as frozen

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FROZEN_SHA = 'cd54e300893a3c657b7b95780523abd37b791032ea691a38b44e389c00129108'
WRAPPER = 'hd_field128_gamma30_d40_r913_midpoint_directed_v25.py'
WRAPPER_SHA = '5c0f188939194cd5f9304644d25a11c8f0b1e97b57a0ce8eddc15f5157751a8b'
SUMMARY_SHA = '3d2a8b5bd5ab38f0417cb9a0d7cd1024374dfd042abbc8c3a3a7368f556cdb40'
PROPOSAL_SHA = '00516791ed2f1bcc112aaadbc97d9e8b3c040ebb58b17cc891cf768410b6ed30'
require, same, fraction = frozen.require, frozen.same, frozen.fraction


def source_pins():
    require(sha256(Path(frozen.__file__).read_bytes()).hexdigest() == FROZEN_SHA,
            'Frozen midpoint auditor changed')
    require(sha256((HERE/WRAPPER).read_bytes()).hexdigest() == WRAPPER_SHA,
            'Reviewed portable producer changed')
    return dict(sorted({**frozen.source_pins(), WRAPPER: WRAPPER_SHA,
                        NAME: sha256((HERE/NAME).read_bytes()).hexdigest()}.items()))


def portable_path(reference):
    """Portable witness paths are theory-relative, not receipt-relative."""
    require(isinstance(reference, dict) and isinstance(reference.get('path'), str)
            and isinstance(reference.get('sha256'), str), 'Missing portable byte link')
    name = reference['path']
    require('\\' not in name and not PurePosixPath(name).is_absolute()
            and not PureWindowsPath(name).is_absolute() and ':' not in name
            and all(part not in ('', '.', '..') for part in name.split('/')),
            'Expected explicit relative portable path')
    path = (HERE/name).resolve()
    path.relative_to(HERE.resolve())
    return path


def public_metadata(record):
    """The frozen public receipt checks, with strict original913 admission."""
    same(record, dict(frozen.previous.EXPECTED, status='complete', failure=None,
        complete_selected_weight=True, covers_all_supports_at_weight=True,
        covers_all_even_odd_splits=True, repetitions=256, paired_region_power=128,
        common_potential_used=False, source_pins_unchanged_at_exit=True,
        source_sha256=frozen.producer.source_pins(), coverage_driver_sha256=frozen.producer.source_pins(),
        bank_source_sha256=frozen.producer.provider.source_pins()), 'complete gamma30/d40 midpoint receipt')
    require(type(record['active_groups']) is int and record['active_groups'] == 913,
            'Expected original913')
    witness = record['witness']
    require(0 < fraction(witness['z']) < 1 and 0 < fraction(witness['alpha']) <= 1,
            'Invalid exact marks')
    require(len(witness['markers']) == 83 and all(len(row) == 2 and all(fraction(x) > 0 for x in row)
            for row in witness['markers']), 'Invalid83 positive marker bank')
    same(witness, dict(state_dimension=30, inner_width=64, gamma_extension_degree=16,
        gamma_base_field_order=256, joint_rank_return_cap=True, all_pair_before_alpha=True,
        precision=256, grid=128, packet_mantissa_bits=80), 'witness scope')
    same(record, dict(expected_even_rows=457, expected_profile_pairs=208849,
        completed_even_rows=457, completed_even_full_range=[0, 456], pending_even_full_range=None,
        profile_pairs_evaluated=208849, requested_failure_bits=54,
        meets_requested_selected_weight_target=True, required_normal_process_exit=True), 'full913 scope')
    require(fraction(record['log2_upper_rational']) <= -54, 'Selected target not met')
    require(record['count_provenance']['maximum_refined_singletons'] == 64,
            'Unexpected near-face cutoff')


def expected_context(record, bank):
    result = {key: value for key, value in bank.items()
              if key not in ('witness_source', 'runner_source_sha256', 'runtime', 'fresh_bank_rebuilt')}
    result.update(source_sha256=frozen.producer.source_pins(),
        coverage_driver_sha256=frozen.producer.source_pins(), witness_source=record['witness_source'],
        active_groups=913, common_potential_used=False,
        transfer_method='frozen positive Arb row-zero sum of (R_even R_odd)^128',
        regional_normalization='frozen RegionalCache, entrywise minimum over all bank markers',
        repetitions=256, paired_region_power=128, requested_failure_bits=54,
        expected_even_rows=457, expected_profile_pairs=208849, max_seconds=record['max_seconds'],
        old_numerical_bounds_reused=False, directed=True, whole_code_certificate=False,
        count_provenance=record['count_provenance'])
    return result


def source_maps(value, pins):
    def merge(mapping):
        require(isinstance(mapping, dict), 'Source map must be an object')
        for name, digest in mapping.items():
            if isinstance(digest, dict):
                merge(digest)
            else:
                require(isinstance(name, str) and isinstance(digest, str)
                        and len(digest) == 64 and all(c in '0123456789abcdef' for c in digest),
                        'Invalid source pin')
                require(name not in pins or pins[name] == digest, 'Conflicting source pin')
                require(sha256((HERE/name).read_bytes()).hexdigest() == digest, 'Source changed: '+name)
                pins[name] = digest
    if isinstance(value, dict):
        for key, child in value.items():
            if key == 'source_sha256' or key.endswith('_source_sha256'):
                merge(child)
            else:
                source_maps(child, pins)
    elif isinstance(value, list):
        for child in value:
            source_maps(child, pins)


def audit_summary(path, *, max_seconds=120.):
    require(type(max_seconds) in (int, float) and 0 < max_seconds <= 120,
            'Expected at most120 seconds')
    started, adapter_pins = monotonic(), source_pins()
    loader, execution_pins = frozen.index.Validator(F(2, 5)), {}
    record, summary_pin = loader.load(path, expected=SUMMARY_SHA, role='portable913-summary')
    public_metadata(record)
    directory = Path(summary_pin['path']).parent
    proposal, proposal_pin = loader.load(portable_path(record['witness_source']),
        expected=PROPOSAL_SHA, role='theory-relative-literal-proposal')
    require(record['witness_source']['sha256'] == PROPOSAL_SHA, 'Wrong proposal link')
    same(record['witness'], proposal['witness'], 'exact positive proposal literals')
    same(proposal, dict(active_groups=913, state_dimension=30, target_distance='2/5',
                       construction=record['construction'], marker_count=83), 'proposal scope')
    source = proposal.get('source')
    if source:
        loader.load(portable_path(source), expected=source['sha256'], role='proposal-source-literals')
    bank, unused = loader.load(directory/'bank.json', role='fresh-bank-metadata')
    driver, unused = loader.load(directory/'driver.json', role='execution-driver-metadata')
    context_file, unused = loader.load(directory/'context.json', role='external-metadata-context')
    same(driver, dict(schema='gamma-r913-midpoint-directed-driver-v25',
        driver=dict(path=WRAPPER, sha256=WRAPPER_SHA), witness_source=record['witness_source'],
        active_groups=913, construction=record['construction'], field_order=1 << 128,
        state_dimension=30, target_distance='2/5', marker_count=83,
        numerical_provider='hd_field128_gamma_midpoint_certificate.py',
        expected_even_rows=457, expected_profile_pairs=208849, precision=256,
        coefficient_bits=80, grid=128, requested_failure_bits=54, maximum_refined_singletons=64,
        whole_job_seconds=600, max_table_bytes=64 << 20, max_output_bytes=2 << 20,
        disk_reserve_bytes=20 << 30, memory_address_space_bytes=4 << 30,
        required_normal_process_exit=True, fresh_bank_required=True, fresh_bank_rebuilt=False,
        old_numerical_bounds_reused=False, directed=False, whole_code_certificate=False), 'actual wrapper scope')
    require(driver['driver'] == record['execution_driver'], 'Wrapper identity mismatch')
    require(driver['source_sha256'] == record['execution_source_sha256'] == bank['runner_source_sha256'],
            'Execution source closure mismatch')
    same(bank, dict(witness_source=record['witness_source'], fresh_bank_rebuilt=True,
                    directed=False, whole_code_certificate=False), 'fresh bank scope')
    require(bank['witness'] == record['witness'], 'Bank/witness mismatch')
    context = expected_context(record, bank)
    digest = sha256(json.dumps(context, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    require(record['context_sha256'] == digest, 'Reconstructed context digest mismatch')
    require(context_file == dict(context=context, context_sha256=digest,
        complete_selected_weight=False, covers_all_supports_at_weight=False, whole_code_certificate=False,
        scope='metadata only; no numerical bound'), 'External context differs from frozen producer metadata')
    require(len(record['rows']) == 457, 'Missing inline rows')
    for position, row in enumerate(record['rows']):
        require(monotonic()-started < max_seconds, 'Saved-data audit budget exhausted')
        require(row['even_full'] == position, 'Inline row order mismatch')
        checkpoint, unused = loader.load(directory/f'row_{position:04d}.json', role='external-even-row')
        require(checkpoint == dict(row, context_sha256=digest, active_groups=913, directed=True,
            whole_code_certificate=False,
            scope='completed even-profile row only; requires midpoint-count metadata context'),
            'External row differs from complete inline row')
    require(len(list(directory.glob('row_*.json'))) == 457, 'Unexpected external row files')
    for value in (record, proposal, bank, driver, context):
        source_maps(value, execution_pins)
    require(monotonic()-started < max_seconds, 'Saved-data audit budget exhausted')
    # Exact frozen arithmetic core: all208849 independent integer count checks,
    # saved row enclosures, a fresh prefix, and the384-bit positive support sum.
    report = frozen.audit_record(record)
    require(monotonic()-started < max_seconds, 'Saved-data audit budget exhausted')
    report.update(receipt=summary_pin, proposal=proposal_pin,
                  requested_failure_bits=54, meets_requested_target=True)
    union = frozen.index.union_ceiling([fraction(report['saved_log2_upper_rational'])], 384)
    require(source_pins() == adapter_pins, 'Adapter closure changed')
    frozen.producer.check_source_pins(execution_pins)
    frozen.producer.check_source_pins(adapter_pins)
    frozen.previous.band._check_links(tuple(loader.files.values()))
    return dict(schema='gamma-midpoint-selected-receipts-audit-v1', audit_result='PASS',
        construction=record['construction'], state_dimension=30, target_distance='2/5',
        selected_weights=[913], weight_audits=[report], precision_bits=384,
        selected_events_union_log2_upper_rational=str(union), source_sha256=frozen.source_pins(),
        checked_file_pins=[summary_pin, proposal_pin], coefficient_replay=False,
        whole_code_certificate=False, coverage_index_modified=False,
        external_checkpoint_context_reconstructed=True, external_rows_checked=457, context_sha256=digest,
        external_input_sha256={Path(pin['path']).relative_to(HERE).as_posix(): pin['sha256']
                               for pin in loader.files.values()},
        external_pin_path_base='Theory directory containing the frozen auditor',
        execution_driver=driver['driver'], execution_source_sha256=dict(sorted(execution_pins.items())),
        audit_adapter_schema='gamma-r913-portable-midpoint-audit-adapter-v1',
        adapter=dict(path=NAME, sha256=adapter_pins[NAME]), adapter_source_sha256=adapter_pins,
        portable_proposal_links_resolved_from_theory=True,
        wrapper_limits_and_exact_proposal_checked=True, external_rows_equal_complete_inline_rows=True,
        source_and_input_bytes_checked_at_exit=True,
        normal_remote_exit_evidence='Parent-reported normal exit; JSON alone does not establish process exit status',
        audit_elapsed_seconds=monotonic()-started,
        audit_limit='Frozen audit_record independently compares every midpoint/final count and recomputes saved-row/prefix384-bit arithmetic. The adapter repeats frozen public metadata checks, resolves portable proposal links from theory, reconstructs the external metadata context and checks every row and execution source byte. Original receipts are unchanged. No coefficient, normalized regional or transfer replay; midpoint theorem and producer arithmetic are inherited from reviewed sources.')
