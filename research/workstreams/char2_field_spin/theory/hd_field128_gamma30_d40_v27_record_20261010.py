"""Publish only two explicitly audited bands through the frozen v27 adapter."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v27 as selected

HERE = Path(__file__).resolve().parent
HELPER_SHA = 'db4d2621547fde7f33ae2aa88a3b360304e12c14778df46751eea59b5ded2c47'
PARENT_SHA = '45d2854157e47e23145f76df445bc8cae4f318ce828003e8322eb36307ea1f20'
RESERVE, WRITE_LIMIT = 20 << 30, 1 << 20
BANDS = (
    ('hd_field128_gamma30_d40_weights209_215_portable_v27_20261010/independent_typed_audit.json',
     '5538711e1cefe9cb30ca7ef6bff0b9534fe368f2da777627be8d27cbd72ea949'),
    ('hd_field128_gamma30_d40_weights1160_1191_portable_v27_20261010/independent_typed_audit.json',
     'a90a2e972389d2e6ce90f88512cdbb014389c874d21d168416349726ff0ba26c'),
)


def reserve_guard(planned_bytes):
    if shutil.disk_usage(HERE).free-planned_bytes < RESERVE:
        raise OSError('Bookkeeping stopped: preserve 20 GiB after planned writes')


def main():
    reserve_guard(WRITE_LIMIT)
    self_path = Path(__file__).resolve()
    self_raw = self_path.read_bytes()
    inputs = {self_path: self_raw}
    parent_path = HERE/'hd_field128_gamma30_d40_partial_cover_v26_20261010.json'
    for path, digest in ((parent_path, PARENT_SHA), (Path(selected.__file__).resolve(), HELPER_SHA),
                         *((HERE/name, digest) for name, digest in BANDS)):
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned assembly input changed: '+str(path))
        inputs[path] = raw
    parent = json.loads(inputs[parent_path])
    result = selected.extend_v26(parent_path, [HERE/name for name, unused in BANDS])
    if (result['receipt_count'] != 255 or result['covered_weight_count'] != 7967
            or result['pending_weight_count'] != 225 or result['added_event_count'] != 2
            or result['added_weight_count'] != 39 or result['entries'][:253] != parent['entries']
            or [row['weight_range'] for row in result['entries'][253:]] != [[209,215],[1160,1191]]
            or result['whole_code_certificate'] is not False
            or result['pending_weight_ranges'] != [[165,179],[521,575],[704,757],[914,999],[1145,1159]]):
        raise ValueError('Unexpected selected ranges, parent prefix, or coverage')
    selection = dict(schema='gamma-v27-explicit-typed-band-selection-v1',
        parent=dict(path=parent_path.name, sha256=PARENT_SHA),
        validator=dict(path=Path(selected.__file__).name, sha256=HELPER_SHA),
        recorder=dict(path=self_path.name, sha256=sha256(self_raw).hexdigest()),
        selected_band_audits=[dict(path=name, sha256=digest) for name,digest in BANDS],
        added_events=2, added_weights=39, previous_253_records_and_order_unchanged=True,
        whole_code_certificate=False, coefficient_replay=False, normalized_regional_replay=False,
        transfer_recomputed=False, output_byte_limit_including_recorder=WRITE_LIMIT,
        minimum_free_bytes_after_planned_writes=RESERVE,
        storage_scope='Existing saved-data bookkeeping only; no coefficient jobs, trees, imports or searches')
    outputs = [(HERE/'hd_field128_gamma30_d40_v27_selection_20261010.json', selection),
               (HERE/'hd_field128_gamma30_d40_partial_cover_v27_20261010.json', result)]
    encoded = [(path, (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()) for path,value in outputs]
    total = len(self_raw)+sum(len(raw) for unused,raw in encoded)
    if total > WRITE_LIMIT:
        raise OSError('Recorder plus selection and index exceed 1 MiB')
    reserve_guard(sum(len(raw) for unused,raw in encoded))
    if any(path.exists() for path,unused in encoded):
        raise FileExistsError('Refusing to replace an immutable index or selection')
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input changed during assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    for path,raw in encoded:
        reserve_guard(len(raw))
        with path.open('xb') as stream:
            stream.write(raw)
        if path.read_bytes() != raw:
            raise RuntimeError('Published bytes differ: '+str(path))
        print(json.dumps(dict(path=str(path), bytes=len(raw), sha256=sha256(raw).hexdigest())), flush=True)
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input changed before closing assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    print(json.dumps(dict(total_new_bytes_including_recorder=total,
        **{key:result[key] for key in ('receipt_count','covered_weight_count','pending_weight_count',
            'pending_weight_ranges','covered_events_union_log2_upper_rational',
            'added_events_union_log2_upper_rational','whole_code_certificate')})), flush=True)


if __name__ == '__main__':
    main()
