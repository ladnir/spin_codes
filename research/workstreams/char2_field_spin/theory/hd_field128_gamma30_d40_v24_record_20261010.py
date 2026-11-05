"""Append four audited remote bands using the frozen v23 index format."""
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v23 as selected
import hd_field128_gamma30_d40_v23_record_20261010 as pilot

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = 'c87acc81cdd8ff1e7128f2f739b78726176247c41045caab773f8119a807a939'
ADDITIONS = (
    ('hd_field128_gamma30_d40_weights140_164_portable_v24_20261010',
     '0bb99f6a55cbc486ed1c7f1649acefd5a14010bce0642b5857e73182263e3f4b', [140, 164]),
    ('hd_field128_gamma30_d40_weights180_208_portable_v24_20261010',
     'c4e95b79d607f094d4a42f6922fe3c987b6445de0f26989062bb62cd0e5b52bd', [180, 208]),
    ('hd_field128_gamma30_d40_weights216_240_portable_v24_20261010',
     'c9ea00f2080d521a16373409f4919ec54e6eb3673b32fb8b4bab95b7dd482842', [216, 240]),
    ('hd_field128_gamma30_d40_weights273_303_portable_nodes2048_v24_20261010',
     '59d7e31c32b0f0694a0928a1990245ededacd19e3eb0d9d4a441c99cabeb6a95', [273, 303]),
)


def main():
    previous_path = HERE/'hd_field128_gamma30_d40_partial_cover_v23_20261010.json'
    previous_raw = previous_path.read_bytes()
    if sha256(previous_raw).hexdigest() != PREVIOUS_SHA:
        raise ValueError('Previous checkpoint changed')
    previous = json.loads(previous_raw)
    inputs = {previous_path: previous_raw, Path(__file__).resolve(): Path(__file__).read_bytes()}
    inputs[Path(pilot.__file__).resolve()] = Path(pilot.__file__).read_bytes()
    if sha256(inputs[Path(pilot.__file__).resolve()]).hexdigest() != 'a890a670c4c9707497c41ac6243062be28e184123557840417044746cc0bdcb4':
        raise ValueError('Pilot selection source changed')
    paths, pins = [], []
    for directory, digest, interval in (*pilot.AUDITS, *ADDITIONS):
        path = HERE/directory/'independent_typed_audit.json'
        raw = path.read_bytes()
        value = json.loads(raw)
        if (sha256(raw).hexdigest() != digest or value['audit_result'] != 'PASS'
                or value['requested_weight_range'] != interval):
            raise ValueError('Selected audited band changed')
        inputs[path] = raw
        paths.append(path)
        pins.append(dict(path=str(path), sha256=digest, weight_range=interval))
    result = selected.extend_v22(HERE/'hd_field128_gamma30_d40_partial_cover_v22_20261010.json', paths)
    if (result['receipt_count'] != 248 or result['covered_weight_count'] != 7728
            or result['pending_weight_count'] != 464 or result['added_weight_count'] != 129
            or result['whole_code_certificate'] is not False
            or result['entries'][:244] != previous['entries']):
        raise ValueError('Unexpected coverage or previous-entry mutation')
    selection = dict(schema='gamma-v24-explicit-remote-band-selection-v1',
        preceding_checkpoint=dict(path=str(previous_path), sha256=PREVIOUS_SHA),
        mathematical_parent=result['parent_index'], selected_since_parent=pins,
        additions_this_round=pins[2:], added_events_this_round=4, added_weights_this_round=110,
        index_format=result['schema'], previous_244_entry_prefix_unchanged=True,
        source_sha256={p.name: sha256(raw).hexdigest() for p, raw in inputs.items() if p.suffix == '.py'},
        scope='Completed audited bands only; partial368..383 and FLOAT913 excluded',
        whole_code_certificate=False)
    outputs = [(HERE/'hd_field128_gamma30_d40_partial_cover_v24_20261010.json', result),
               (HERE/'hd_field128_gamma30_d40_v24_selection_20261010.json', selection)]
    encoded = [(p, (json.dumps(v, indent=2, allow_nan=False)+'\n').encode()) for p, v in outputs]
    size = sum(len(raw) for p, raw in encoded)
    if size > 5 << 20 or shutil.disk_usage(HERE).free-size < 1 << 30:
        raise OSError('Bounded bookkeeping budget exceeded')
    if any(p.exists() for p, raw in encoded):
        raise FileExistsError('Refusing to replace an index checkpoint')
    if any(p.read_bytes() != raw for p, raw in inputs.items()):
        raise RuntimeError('Input changed before output')
    for p, raw in encoded:
        if shutil.disk_usage(HERE).free-len(raw) < 1 << 30:
            raise OSError('Bookkeeping reserve exhausted')
        with p.open('xb') as f:
            f.write(raw)
        print(json.dumps(dict(path=str(p), sha256=sha256(raw).hexdigest(), bytes=len(raw))), flush=True)
    print(json.dumps({k: result[k] for k in ('covered_weight_count', 'pending_weight_count',
        'receipt_count', 'pending_weight_ranges', 'covered_events_union_log2_upper_rational')}), flush=True)


if __name__ == '__main__':
    main()
