"""Record the two returned, independently audited portable pilot bands."""
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v23 as selected

HERE = Path(__file__).resolve().parent
AUDITS = (
    ('hd_field128_gamma30_d40_weights65_71_portable_v23_20261010',
     '555c0c4a5335d6cff2a84b97055c832377e4421709681eea9a3daeee73963e5a', [65, 71]),
    ('hd_field128_gamma30_d40_weights89_100_portable_v23_20261010',
     'ed0923ed285ba7fc4578f2219417c82c718183b7966a1cf2d3b868d7798b0a6f', [89, 100]),
)


def main():
    source = Path(__file__).resolve()
    raw_source = source.read_bytes()
    paths, inputs, additions = [], {}, []
    for directory, digest, interval in AUDITS:
        path = HERE / directory / 'independent_typed_audit.json'
        raw = path.read_bytes()
        row = json.loads(raw)
        if (sha256(raw).hexdigest() != digest or row['audit_result'] != 'PASS'
                or row['requested_weight_range'] != interval):
            raise ValueError('Pilot audit identity differs')
        inputs[path], paths = raw, paths + [path]
        additions.append(dict(path=str(path), sha256=digest, weight_range=interval))
    parent = HERE / 'hd_field128_gamma30_d40_partial_cover_v22_20261010.json'
    result = selected.extend_v22(parent, paths)
    if (result['receipt_count'] != 244 or result['covered_weight_count'] != 7618
            or result['pending_weight_count'] != 574 or result['added_weight_count'] != 19
            or result['whole_code_certificate'] is not False):
        raise ValueError('Unexpected pilot coverage')
    selection = dict(schema='gamma-v23-explicit-portable-pilot-selection-v1',
        parent=dict(path=str(parent), sha256=selected.PARENT_SHA256), additions=additions,
        added_event_count=2, added_weight_count=19,
        source_sha256={source.name: sha256(raw_source).hexdigest()},
        scope='Only the two complete returned bands; no unresolved or partial attempt admitted',
        whole_code_certificate=False)
    outputs = (
        (HERE / 'hd_field128_gamma30_d40_v23_selection_20261010.json', selection),
        (HERE / 'hd_field128_gamma30_d40_partial_cover_v23_20261010.json', result),
    )
    encoded = [(path, (json.dumps(value, indent=2, allow_nan=False)+'\n').encode())
               for path, value in outputs]
    size = sum(len(raw) for path, raw in encoded)
    if (size > 5 << 20 or shutil.disk_usage(HERE).free-size < 600 << 20
            or any(path.exists() for path, raw in encoded)):
        raise RuntimeError('Bounded exclusive record preflight failed')
    if source.read_bytes() != raw_source or any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Record input changed')
    for path, raw in encoded:
        if shutil.disk_usage(HERE).free-len(raw) < 600 << 20:
            raise RuntimeError('Record write would breach local reserve')
        with path.open('xb') as stream:
            stream.write(raw)
        print(json.dumps(dict(path=str(path), sha256=sha256(raw).hexdigest())), flush=True)
    print(json.dumps({key: result[key] for key in ('receipt_count', 'covered_weight_count',
        'pending_weight_count', 'covered_events_union_log2_upper_rational',
        'added_events_union_log2_upper_rational', 'pending_weight_ranges')}), flush=True)


if __name__ == '__main__':
    main()
