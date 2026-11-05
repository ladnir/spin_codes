"""Assemble exactly 37 audited additions while preserving the pinned v20 index."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_selected_index_v21 as index


HERE = Path(__file__).resolve().parent
MANIFESTS = (
    ('hd_field128_gamma30_d40_kernel_v21_band_manifest_20261010.json',
     'f68f3cd4b8c7cf2a7ad5c1574cae2e2ab1be875e77ed028d858fd7c6e54ecaf8',
     'events', 26),
    ('hd_field128_gamma30_d40_gaps2257_2431_2848_3007_manifest_20261010.json',
     '32cb1d5dd46462730b5e00a70afd198cb6f65edecd44beb7f3ab309403bd1545',
     'individually_audited_events', 6),
)
ROOT_BANDS = (
    (101, 103, 'r112', '213788c3d9fccf8ac9fa82c10e33fcf67c9d3efe929e9c47e27c1126e4315504'),
    (104, 120, 'r112', '72b9116f51b02383f04f47223c5219cfba22cbd71ef476bd6b55096a9d45e46c'),
    (121, 127, 'r128', 'a2b088813ff20daa2fa8d29af4aeb3bb7cd60d81ddc589d119da24eec955c558'),
    (129, 139, 'r128', 'bacb307a60cdb6b6c35f67b0251216083a2d87b7c544495f6e71b0bcf08d1129'),
    (3904, 3936, 'r3936', '41604929f57a858b359758750b57a65ec64207551f8cb2cfaf0098eb25e87c79'),
)


def load(path, expected=None):
    path = Path(path).resolve()
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if expected is not None and digest != expected:
        raise RuntimeError('Immutable input changed: '+str(path))
    return json.loads(raw), dict(path=str(path), sha256=digest)


def main():
    source = Path(__file__).resolve()
    source_hash = sha256(source.read_bytes()).hexdigest()
    manifests, audits, extra_inputs = [], [], []
    for name, digest, key, count in MANIFESTS:
        manifest, pin = load(HERE/name, digest)
        if len(manifest[key]) != count:
            raise ValueError('Unexpected selected event count')
        manifests.append(pin)
        for event in manifest[key]:
            unused, audit_pin = load(event['independent_audit']['path'],
                                     event['independent_audit']['sha256'])
            audits.append(audit_pin)
            if 'execution_driver_receipt' in event:
                receipt = event['execution_driver_receipt']
                unused, receipt_pin = load(receipt['path'], receipt['sha256'])
                extra_inputs.append(receipt_pin)
    for lo, hi, tag, digest in ROOT_BANDS:
        folder = HERE/f'hd_field128_gamma30_d40_weights{lo}_{hi}_{tag}_band_v21_20261010'
        audit, pin = load(folder/'independent_audit.json', digest)
        if audit['requested_weight_range'] != [lo, hi]:
            raise ValueError('Unexpected root audit interval')
        audits.append(pin)
    if len(audits) != 37 or len({p['path'] for p in audits}) != 37:
        raise ValueError('Selection is not 37 distinct audited events')
    result = index.extend_v20(HERE/'hd_field128_gamma30_d40_partial_cover_v20_20261010.json',
                             [p['path'] for p in audits])
    if (result['added_event_count'] != 37 or result['added_weight_count'] != 853
            or result['covered_weight_count'] != 7447 or result['pending_weight_count'] != 745):
        raise ValueError('Unexpected new coverage')
    for pin in [*manifests, *audits, *extra_inputs]:
        load(pin['path'], pin['sha256'])
    if sha256(source.read_bytes()).hexdigest() != source_hash:
        raise RuntimeError('Assembly source changed')
    selection = dict(schema='gamma-partial-cover-v21-explicit-selection-v1',
        construction=result['construction'], target_distance='2/5',
        manifests=manifests, root_bands=[list(row[:3]) for row in ROOT_BANDS],
        selected_audits=audits, selected_event_count=37, selected_weight_count=853,
        extra_execution_driver_receipts=extra_inputs,
        assembly_source=dict(path=str(source), sha256=source_hash),
        partial_attempts_not_selected=True, no_weight_interpolation=True,
        coefficient_replay=False, whole_code_certificate=False)
    selection_pin = index.frozen.producer.profiles.write_new(
        HERE/'hd_field128_gamma30_d40_v21_selection_20261010.json', selection)
    result_pin = index.frozen.producer.profiles.write_new(
        HERE/'hd_field128_gamma30_d40_partial_cover_v21_20261010.json', result)
    print(json.dumps(dict(selection=selection_pin, index=result_pin,
        covered=result['covered_weight_count'], pending=result['pending_weight_count'],
        events=result['receipt_count'], union=result['covered_events_union_log2_upper_rational'],
        new_union=result['added_events_union_log2_upper_rational'],
        pending_ranges=result['pending_weight_ranges'])), flush=True)


if __name__ == '__main__':
    main()
