"""Append the audited368..383 band using the frozen v23 index validator."""
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v23 as selected
import hd_field128_gamma30_d40_v23_record_20261010 as pilot
import hd_field128_gamma30_d40_v24_record_20261010 as prior

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = '703252867f063c8bff2b76f83c582dbb8675577e6a2f56d7e07093b6430d5c45'
SOURCE_PINS = {
    'hd_field128_gamma_selected_index_v23.py':
        '5543d4984d341204cc8de168a03a18cc1d7dbd700e3890debe097de7305808c8',
    'hd_field128_gamma30_d40_v23_record_20261010.py':
        'a890a670c4c9707497c41ac6243062be28e184123557840417044746cc0bdcb4',
    'hd_field128_gamma30_d40_v24_record_20261010.py':
        'b59e5b688da796272b6f19f001b45ecf8904c2705640f8032a3eba2e5a0431db',
}
ADDITIONS = (
    ('hd_field128_gamma30_d40_weights368_383_portable_nodes2048_v25_20261010',
     '52fb710b3834b0f6b3b53e095bceb0fefd857afb1d86a8e01f02715ca8262d56', [368, 383]),
)
MAX_OUTPUT_BYTES, RESERVE = 2 << 20, 20 << 30


def main():
    previous_path = HERE/'hd_field128_gamma30_d40_partial_cover_v24_20261010.json'
    previous_raw = previous_path.read_bytes()
    if sha256(previous_raw).hexdigest() != PREVIOUS_SHA:
        raise ValueError('Previous checkpoint changed')
    previous = json.loads(previous_raw)
    if (previous['receipt_count'] != 248 or len(previous['entries']) != 248
            or previous['covered_weight_count'] != 7728
            or previous['pending_weight_count'] != 464
            or previous['whole_code_certificate'] is not False):
        raise ValueError('Unexpected preceding checkpoint scope')
    inputs = {previous_path: previous_raw, Path(__file__).resolve(): Path(__file__).read_bytes()}
    for name, expected in SOURCE_PINS.items():
        path = HERE/name
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != expected:
            raise ValueError('Frozen selection source changed: '+name)
        inputs[path] = raw
    paths, pins = [], []
    for directory, digest, interval in (*pilot.AUDITS, *prior.ADDITIONS, *ADDITIONS):
        path = HERE/directory/'independent_typed_audit.json'
        raw = path.read_bytes()
        value = json.loads(raw)
        if (sha256(raw).hexdigest() != digest or value['audit_result'] != 'PASS'
                or value['requested_weight_range'] != interval):
            raise ValueError('Selected audited band changed')
        inputs[path] = raw
        paths.append(path)
        pins.append(dict(path=str(path), sha256=digest, weight_range=interval))
    if len(paths) != 7:
        raise ValueError('Expected exactly seven audited events since v22')
    result = selected.extend_v22(HERE/'hd_field128_gamma30_d40_partial_cover_v22_20261010.json', paths)
    if (result['receipt_count'] != 249 or result['covered_weight_count'] != 7744
            or result['pending_weight_count'] != 448 or result['added_weight_count'] != 145
            or result['whole_code_certificate'] is not False
            or result['entries'][:248] != previous['entries']):
        raise ValueError('Unexpected coverage or previous-entry mutation')
    selection = dict(schema='gamma-v25-explicit-remote-band-selection-v1',
        preceding_checkpoint=dict(path=str(previous_path), sha256=PREVIOUS_SHA),
        mathematical_parent=result['parent_index'], selected_since_parent=pins,
        additions_this_round=pins[6:], added_events_this_round=1, added_weights_this_round=16,
        index_format=result['schema'], previous_248_entry_prefix_unchanged=True,
        source_sha256={p.name: sha256(raw).hexdigest() for p, raw in inputs.items() if p.suffix == '.py'},
        scope='Completed audited368..383 only this round; earlier partials and913 excluded',
        whole_code_certificate=False)
    outputs = [(HERE/'hd_field128_gamma30_d40_partial_cover_v25_20261010.json', result),
               (HERE/'hd_field128_gamma30_d40_v25_selection_20261010.json', selection)]
    encoded = [(p, (json.dumps(v, indent=2, allow_nan=False)+'\n').encode()) for p, v in outputs]
    size = sum(len(raw) for p, raw in encoded)
    if size > MAX_OUTPUT_BYTES or shutil.disk_usage(HERE).free-size < RESERVE:
        raise OSError('Bounded bookkeeping budget or20GiB reserve exceeded')
    if any(p.exists() for p, raw in encoded):
        raise FileExistsError('Refusing to replace an index checkpoint')
    if any(p.read_bytes() != raw for p, raw in inputs.items()):
        raise RuntimeError('Input changed before output')
    for p, raw in encoded:
        if shutil.disk_usage(HERE).free-len(raw) < RESERVE:
            raise OSError('Bookkeeping reserve exhausted')
        with p.open('xb') as f:
            f.write(raw)
        print(json.dumps(dict(path=str(p), sha256=sha256(raw).hexdigest(), bytes=len(raw))), flush=True)
    print(json.dumps({k: result[k] for k in ('covered_weight_count', 'pending_weight_count',
        'receipt_count', 'pending_weight_ranges', 'covered_events_union_log2_upper_rational')}), flush=True)


if __name__ == '__main__':
    main()
