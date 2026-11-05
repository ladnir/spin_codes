"""Publish four compact remote-evidence pairs; never import or replay their trees."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v28 as selected

HERE = Path(__file__).resolve().parent
HELPER_SHA = 'ff1d62fbeee9b232319cea920a26a7485bacb735ce0aa9ae59b006e2dc1246d2'
TEST_SHA = 'e8951b781eafe71175f8fbeda51b0f067f48e44314a82eec4e34020c28365e93'
PARENT_SHA = 'bb5873cc96c861d3c9f08806bdbae63c5778fd9d3d57999a2631f529dd5e3a5f'
CHECK_SHA = '2bd521821963ef03eeeec5ad55ff0945b0a6ccb90ca616e0842e4dda8e80dead'
RESERVE, WRITE_LIMIT, AUDIT_ALLOWANCE = 20 << 30, 1 << 20, 16 << 10
PAIRS = (
    ('165_179', 'ff3779a7fe84937c8e5cbcbd512a8b426c30474a1c62366b09d0e1050293c7c6',
     'df41ed601d582909a60e7e3ba06cde28274a9a2f93ab97067f1404a8d21e74c6'),
    ('704_757', '87abfead338b9462897071807fb5a73f6ec8b52bffde20c4404bfcc5da883627',
     '5cfc33a8163ae4e5c74e0936d52cf1bfa9bf244b1d57b64cc3f6acf99ac412b3'),
    ('968_999', 'a75512816ac9d75e34da2d77eb5b4083d6a98b218f4c09fd0737ec8cfd8638f0',
     '680a6575245d2dbbccfe427e40781d13664147bba370ecffb3c5e126e01d0a59'),
    ('1150_1159', '10e5f38ca816f685d544268972a2d31df3b611dbb99649dc4a3a75175aa4fc52',
     '551083d7fee6d4856c0ef93fa949a76ada3ec04d48114f97d34ea85950fb6c82'),
)


def event_specs():
    return [dict(audit=dict(path=f'remote_receipts_v28/audit_{name}.json', sha256=audit),
                 archive_receipt=dict(path=f'remote_receipts_v28/archive_{name}.json', sha256=archive))
            for name, audit, archive in PAIRS]


def reserve_guard(planned_bytes):
    if shutil.disk_usage(HERE).free-planned_bytes < RESERVE:
        raise OSError('Bookkeeping stopped: preserve 20 GiB after planned writes')


def main():
    reserve_guard(WRITE_LIMIT)
    self_path = Path(__file__).resolve()
    self_raw = self_path.read_bytes()
    inputs = {self_path: self_raw}
    parent_path = HERE/'hd_field128_gamma30_d40_partial_cover_v27_20261010.json'
    check_path = HERE/'hd_field128_gamma30_d40_partial_cover_v27_20261010_independent_check.json'
    specs = event_specs()
    explicit = [(parent_path, PARENT_SHA), (check_path, CHECK_SHA),
                (Path(selected.__file__).resolve(), HELPER_SHA),
                (HERE/'hd_field128_gamma_selected_index_v28_tests.py', TEST_SHA)]
    explicit += [(HERE/pin['path'], pin['sha256']) for pair in specs for pin in pair.values()]
    for path, digest in explicit:
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned assembly input changed: '+str(path))
        inputs[path] = raw
    parent = json.loads(inputs[parent_path])
    result = selected.extend_v27(parent_path, check_path, specs)
    if (result['receipt_count'] != 259 or result['covered_weight_count'] != 8078
            or result['pending_weight_count'] != 114 or result['added_event_count'] != 4
            or result['added_weight_count'] != 111 or result['entries'][:255] != parent['entries']
            or [row['weight_range'] for row in result['entries'][255:]] != [[165,179],[704,757],[968,999],[1150,1159]]
            or result['pending_weight_ranges'] != [[521,575],[914,967],[1145,1149]]
            or result['whole_code_certificate'] is not False):
        raise ValueError('Unexpected selected ranges, parent prefix, or coverage')
    selection = dict(schema='gamma-v28-explicit-remote-evidence-selection-v1',
        parent=dict(path=parent_path.name, sha256=PARENT_SHA),
        parent_independent_check=dict(path=check_path.name, sha256=CHECK_SHA),
        validator=dict(path=Path(selected.__file__).name, sha256=HELPER_SHA),
        validator_tests=dict(path='hd_field128_gamma_selected_index_v28_tests.py', sha256=TEST_SHA),
        recorder=dict(path=self_path.name, sha256=sha256(self_raw).hexdigest()),
        event_specs=specs, added_events=4, added_weights=111,
        previous_255_records_and_order_unchanged=True,
        mathematical_audits_inherited_from_remote=True, archive_verification_inherited_from_remote=True,
        historical_origin_verified=False, local_tree_replay=False, local_archive_revalidation=False,
        whole_code_certificate=False, coefficient_replay=False, normalized_regional_replay=False,
        transfer_recomputed=False, output_byte_limit_including_sources_and_independent_audit=WRITE_LIMIT,
        independent_audit_byte_allowance=AUDIT_ALLOWANCE, minimum_free_bytes_after_planned_writes=RESERVE,
        storage_scope='Compact local records only; remote tree/archive verification inherited, not repeated locally')
    outputs = [(HERE/'hd_field128_gamma30_d40_v28_selection_20261010.json', selection),
               (HERE/'hd_field128_gamma30_d40_partial_cover_v28_20261010.json', result)]
    encoded = [(path,(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()) for path,value in outputs]
    sources = [self_path,Path(selected.__file__).resolve(),HERE/'hd_field128_gamma_selected_index_v28_tests.py']
    total = sum(len(inputs[path]) for path in sources)+sum(len(raw) for unused,raw in encoded)+AUDIT_ALLOWANCE
    if total > WRITE_LIMIT:
        raise OSError('Sources, selection, index and reserved audit exceed one MiB')
    reserve_guard(sum(len(raw) for unused,raw in encoded)+AUDIT_ALLOWANCE)
    if any(path.exists() for path,unused in encoded):
        raise FileExistsError('Refusing to replace an immutable index or selection')
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input changed during assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    for path,raw in encoded:
        reserve_guard(len(raw)+AUDIT_ALLOWANCE)
        with path.open('xb') as stream:
            stream.write(raw)
        if path.read_bytes() != raw:
            raise RuntimeError('Published bytes differ: '+str(path))
        print(json.dumps(dict(path=str(path),bytes=len(raw),sha256=sha256(raw).hexdigest())),flush=True)
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input changed before closing assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    print(json.dumps(dict(planned_bytes_including_sources_and_reserved_audit=total,
        **{key:result[key] for key in ('receipt_count','covered_weight_count','pending_weight_count',
            'pending_weight_ranges','covered_events_union_log2_upper_rational',
            'added_events_union_log2_upper_rational','whole_code_certificate')})),flush=True)


if __name__ == '__main__':
    main()
