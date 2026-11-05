"""Publish an explicitly selected, audited extension; no new numerical caps."""
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v26 as selected

HERE = Path(__file__).resolve().parent
HELPER_SHA = 'fdaac0de6128eb17007a811611dce52aa2c8c356ab4037c9e99a6c63ec2bb0f2'
PARENT_SHA = 'c92b67290ff9b1a729496598e071a7599a3c76ce2486e6441023ab3dee59aea0'
POINT = ('hd_field128_gamma30_d40_weights913_913_midpoint_v25_20261010/independent_audit.json',
         '6af8461d0ebd10a2ca7432f13b809c2e28ab72cab4ce69dc0a4d5daac526c6cb')
BANDS = (
    ('hd_field128_gamma30_d40_weights385_439_portable_v26_20261010/independent_typed_audit.json',
     'e62c3ec0a559cb095f9b59def1ff9040c697635adf4236dd08adf39e3579e06f'),
    ('hd_field128_gamma30_d40_weights576_639_portable_v26_20261010/independent_typed_audit.json',
     '7f3a8ad69d6e9446bc064276a1ba936243e75d19203d0b6636ce088d75287a8c'),
    ('hd_field128_gamma30_d40_weights640_703_portable_v26_20261010/independent_typed_audit.json',
     '0322fc7ee4c3616429cadf463897203e9e9e642a9386b347b0aaceedd8d64565'),
)


def main():
    parent_path = HERE/'hd_field128_gamma30_d40_partial_cover_v25_20261010.json'
    inputs = {}
    for path, expected in ((parent_path, PARENT_SHA), (Path(selected.__file__), HELPER_SHA),
                           *((HERE/name, digest) for name, digest in (POINT, *BANDS))):
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != expected:
            raise ValueError('Pinned input changed: '+str(path))
        inputs[path] = raw
    parent = json.loads(inputs[parent_path])
    result = selected.extend_v25(parent_path, HERE/POINT[0], [HERE/name for name, unused in BANDS])
    if (result['receipt_count'] != 253 or result['covered_weight_count'] != 7928
            or result['pending_weight_count'] != 264 or result['added_weight_count'] != 184
            or result['entries'][:249] != parent['entries']
            or result['whole_code_certificate'] is not False
            or result['covered_events_union_log2_upper_rational'] != '-6846934658051/125000000000'
            or result['pending_weight_ranges'] != [[165,179],[209,215],[521,575],
                                                    [704,757],[914,999],[1145,1191]]):
        raise ValueError('Unexpected scope, parent mutation, or union')
    selection = dict(schema='gamma-v26-explicit-mixed-event-selection-v1',
        parent=dict(path=parent_path.name, sha256=PARENT_SHA),
        validator=dict(path=Path(selected.__file__).name, sha256=HELPER_SHA),
        recorder=dict(path=Path(__file__).name, sha256=sha256(Path(__file__).read_bytes()).hexdigest()),
        selected_midpoint_audit=dict(path=POINT[0], sha256=POINT[1]),
        selected_band_audits=[dict(path=name, sha256=digest) for name,digest in BANDS],
        added_events=4, added_weights=184, previous_249_records_and_order_unchanged=True,
        whole_code_certificate=False, coefficient_replay=False,
        output_byte_limit=1 << 20,
        storage_scope='Bounded maintenance of existing results under AGENTS.md; no search, coefficient generation, or build')
    outputs = [(HERE/'hd_field128_gamma30_d40_partial_cover_v26_20261010.json', result),
               (HERE/'hd_field128_gamma30_d40_v26_selection_20261010.json', selection)]
    encoded = [(path, (json.dumps(value, indent=2, allow_nan=False)+'\n').encode())
               for path,value in outputs]
    total = sum(len(raw) for unused,raw in encoded)
    if total > 1 << 20 or shutil.disk_usage(HERE).free-total < 16 << 20:
        raise OSError('Bounded bookkeeping write budget or available space exhausted')
    if any(path.exists() for path,unused in encoded):
        raise FileExistsError('Refusing to replace an index checkpoint')
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input bytes changed during assembly')
    for path,raw in encoded:
        with path.open('xb') as stream:
            stream.write(raw)
        print(json.dumps(dict(path=str(path), sha256=sha256(raw).hexdigest(), bytes=len(raw))), flush=True)
    print(json.dumps({key:result[key] for key in ('receipt_count','covered_weight_count',
        'pending_weight_count','pending_weight_ranges','covered_events_union_log2_upper_rational')}), flush=True)


if __name__ == '__main__':
    main()
