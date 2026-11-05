"""Publish nine pinned remote audits; a separate completion audit remains required."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re
import shutil

import hd_field128_gamma_selected_index_v31 as selected

HERE = Path(__file__).resolve().parent
HELPER_SHA = 'a98322ce10ae59a98bd58a8ba6a634058402dd40111fab5d7ccf5bf901dc7328'
TEST_SHA = '7efa1e843a513a7d916755341b9029989b36fd939a2ecdc7bf192f4766824604'
PARENT_SHA = '875c572f4ae84b5bc139a4c21be9b0b4e18aed51ab1f288c45c125a260aab88f'
CHECK_SHA = '22357f681f345f94f13d72394df478984f8756d5b68fbaf0db917d0b6e9b604a'
RESERVE, OUTPUT_LIMIT, RECORDER_LIMIT = 20 << 30, 512 << 10, 12 << 10
EXPECTED_RANGES = ((521, 547), (549, 575), (917, 917), (918, 918),
                   (919, 919), (1146, 1146), (1147, 1147), (1148, 1148), (1149, 1149))


def reserve_guard(planned_bytes):
    if shutil.disk_usage(HERE).free - planned_bytes < RESERVE:
        raise OSError('Bookkeeping stopped: preserve 20 GiB after planned writes')


def local_path(name):
    if not isinstance(name, str) or '\\' in name:
        raise ValueError('Expected a theory-relative POSIX path')
    path = PurePosixPath(name)
    if path.is_absolute() or not path.parts or any(p in ('.', '..') or ':' in p for p in path.parts):
        raise ValueError('Expected a non-escaping relative path')
    resolved = (HERE / path).resolve()
    if not resolved.is_relative_to(HERE) or resolved == HERE:
        raise ValueError('Input escapes theory directory')
    return resolved


def digest(value):
    if not isinstance(value, str) or re.fullmatch('[0-9a-f]{64}', value) is None:
        raise ValueError('Expected an explicit lowercase SHA256')
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection-source', required=True,
                        help='Pinned compact JSON containing exactly nine event_specs')
    parser.add_argument('--selection-sha256', required=True)
    args = parser.parse_args(argv)
    reserve_guard(OUTPUT_LIMIT + RECORDER_LIMIT)
    self_path = Path(__file__).resolve()
    self_raw = self_path.read_bytes()
    if len(self_raw) > RECORDER_LIMIT:
        raise OSError('Recorder source exceeds 12 KiB')
    inputs = {self_path: self_raw}

    def read_pinned(path, expected):
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest(expected):
            raise ValueError('Pinned assembly input changed: ' + str(path))
        inputs[path] = raw
        return raw

    proposal_path = local_path(args.selection_source)
    proposal_raw = read_pinned(proposal_path, args.selection_sha256)
    proposal = json.loads(proposal_raw)
    specs = proposal.get('event_specs') if isinstance(proposal, dict) else None
    if not isinstance(specs, list) or len(specs) != 9:
        raise ValueError('Expected exactly nine explicit audit/archive pairs')
    input_names = set()
    for spec in specs:
        if not isinstance(spec, dict) or set(spec) != {'audit', 'archive_receipt'}:
            raise ValueError('Expected audit/archive_receipt pair')
        for pin in spec.values():
            if not isinstance(pin, dict) or set(pin) != {'path', 'sha256'}:
                raise ValueError('Expected explicit path/SHA256')
            path = local_path(pin['path'])
            if path in input_names:
                raise ValueError('Duplicate compact receipt path')
            input_names.add(path)
            read_pinned(path, pin['sha256'])

    parent_path = HERE / 'hd_field128_gamma30_d40_partial_cover_v30_20261010.json'
    check_path = HERE / 'hd_field128_gamma30_d40_partial_cover_v30_20261010_independent_check.json'
    for path, expected in ((parent_path, PARENT_SHA), (check_path, CHECK_SHA),
                          (Path(selected.__file__).resolve(), HELPER_SHA),
                          (HERE / 'hd_field128_gamma_selected_index_v31_tests.py', TEST_SHA)):
        read_pinned(path, expected)
    parent = json.loads(inputs[parent_path])
    result = selected.extend_v30(parent_path, check_path, specs)
    added = result['entries'][265:]
    actual_ranges = [tuple(row['weight_range']) for row in added]
    if (result['receipt_count'] != 274 or result['covered_weight_count'] != 8192
            or result['pending_weight_count'] != 0 or result['added_event_count'] != 9
            or result['added_weight_count'] != 61
            or selected.canonical(result['entries'][:265]) != selected.canonical(parent['entries'])
            or len(actual_ranges) != len(set(actual_ranges))
            or set(actual_ranges) != set(EXPECTED_RANGES)
            or result['covered_weight_ranges'] != [[1, 8192]]
            or result['pending_weight_ranges'] != []
            or result['covered_events_failure_bits'] != 40
            or result['aggregation_precision_bits'] != 384
            or result['meets_covered_events_target'] is not True
            or F(result['covered_events_union_log2_upper_rational']) > -40
            or result['whole_code_certificate'] is not False
            or result['completion_audit_performed'] is not False
            or len(result['latest_audited_multiweight_events']) != 2
            or len(result['latest_explicit_selected_events']) != 7):
        raise ValueError('Unexpected events, parent prefix, coverage, union or completion scope')
    selection = dict(schema='gamma-v31-explicit-remote-evidence-selection-v1',
        input_selection=dict(path=proposal_path.relative_to(HERE).as_posix(),
                             sha256=args.selection_sha256),
        parent=dict(path=parent_path.name, sha256=PARENT_SHA),
        parent_independent_check=dict(path=check_path.name, sha256=CHECK_SHA),
        validator=dict(path=Path(selected.__file__).name, sha256=HELPER_SHA),
        validator_tests=dict(path='hd_field128_gamma_selected_index_v31_tests.py', sha256=TEST_SHA),
        recorder=dict(path=self_path.name, sha256=sha256(self_raw).hexdigest()),
        event_specs=specs, added_events=9, added_weights=61,
        previous_265_records_and_order_unchanged=True,
        added_weight_ranges=[list(pair) for pair in actual_ranges],
        actual_ceiling_union_failure_bits=40, individual_requested_targets_preserved=True,
        composite_bands_included=True, composite_band_count=2, selected_singleton_count=7,
        nonoverlapping_weight_ranges=True, mathematical_audits_inherited_from_remote=True,
        archive_verification_inherited_from_remote=True, historical_origin_verified=False,
        local_tree_replay=False, local_archive_revalidation=False, coefficient_replay=False,
        normalized_regional_replay=False, transfer_recomputed=False,
        all_weight_classes_indexed=True, completion_audit_performed=False,
        whole_code_certificate=False, combined_selection_index_byte_limit=OUTPUT_LIMIT,
        recorder_source_byte_limit=RECORDER_LIMIT, minimum_free_bytes_after_planned_writes=RESERVE,
        storage_scope='Compact records only; remote audits and durable archive verification are inherited',
        completion_scope='All weights indexed; separate completion audit required, no whole-code claim')
    outputs = [(HERE / 'hd_field128_gamma30_d40_v31_selection_20261010.json', selection),
               (HERE / 'hd_field128_gamma30_d40_partial_cover_v31_20261010.json', result)]
    encoded = [(path, (json.dumps(value, indent=2, allow_nan=False) + '\n').encode())
               for path, value in outputs]
    total = sum(len(raw) for unused, raw in encoded)
    if total > OUTPUT_LIMIT:
        raise OSError('Selection and index exceed 512 KiB')
    reserve_guard(total)
    if any(path.exists() for path, unused in encoded):
        raise FileExistsError('Refusing to replace immutable selection or index')
    if any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Input changed during assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    remaining = total
    for path, raw in encoded:
        reserve_guard(remaining)
        with path.open('xb') as stream:
            stream.write(raw)
        remaining -= len(raw)
        if path.read_bytes() != raw:
            raise RuntimeError('Published bytes differ: ' + str(path))
        print(json.dumps(dict(path=str(path), bytes=len(raw), sha256=sha256(raw).hexdigest())), flush=True)
    if any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Input changed before closing assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    print(json.dumps(dict(selection_index_bytes=total, recorder_bytes=len(self_raw),
        **{key: result[key] for key in ('receipt_count', 'covered_weight_count', 'pending_weight_count',
            'pending_weight_ranges', 'covered_events_union_log2_upper_rational',
            'added_events_union_log2_upper_rational', 'whole_code_certificate',
            'completion_audit_performed')})), flush=True)


if __name__ == '__main__':
    main()
