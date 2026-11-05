"""Assemble exactly15 independently audited additions; retain frozen v21."""
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v22_typed as index

HERE = Path(__file__).resolve().parent
ROOT_BANDS = [(72, 79, 'r80'), (80, 88, 'r80')]
ROOT_BANDS += [(lo, lo+7, 'r320') for lo in range(304, 368, 8)]


def load(path, expected=None):
    path = Path(path).resolve()
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if expected is not None and digest != expected:
        raise RuntimeError('Immutable input changed: '+str(path))
    return json.loads(raw), dict(path=str(path), sha256=digest)


def write_new(path, data):
    # Final bookkeeping is bounded separately from numerical checkpoint runs.
    # Existing numerical producers retain their750MiB reserve unchanged.
    if shutil.disk_usage(HERE).free < 600 << 20:
        raise RuntimeError('Final bookkeeping reserve below600MiB')
    if len(json.dumps(data).encode()) > 5 << 20:
        raise RuntimeError('Unexpected bookkeeping output above5MiB')
    return index.producer.profiles.write_new(path, data)


def main():
    own = Path(__file__).resolve()
    own_hash = sha256(own.read_bytes()).hexdigest()
    manifest, manifest_pin = load(HERE/'hd_field128_gamma30_d40_kernel_v22_band_manifest_20261010.json')
    if len(manifest['events']) != 5:
        raise ValueError('Expected five explicitly selected kernel events')
    audits = []
    for event in manifest['events']:
        row, pin = load(event['independent_audit']['path'], event['independent_audit']['sha256'])
        audits.append(pin)
    actual_ranges = sorted(load(p['path'], p['sha256'])[0]['requested_weight_range'] for p in audits)
    if actual_ranges != [[440, 447], [448, 463], [464, 479], [480, 495], [496, 510]]:
        raise ValueError('Kernel selection is not exactly440..510 in its five completed bands')
    for lo, hi, tag in ROOT_BANDS:
        directory = HERE/f'hd_field128_gamma30_d40_weights{lo}_{hi}_{tag}_v22_20261010'
        row, pin = load(directory/'independent_typed_audit.json')
        if row['requested_weight_range'] != [lo, hi]:
            raise ValueError('Unexpected explicit root audit interval')
        audits.append(pin)
    if len(audits) != 15 or len({p['path'] for p in audits}) != 15:
        raise ValueError('Expected fifteen distinct audits')
    result = index.extend_v21(HERE/'hd_field128_gamma30_d40_partial_cover_v21_20261010.json',
                             [p['path'] for p in audits])
    if (result['added_event_count'] != 15 or result['added_weight_count'] != 152
            or result['covered_weight_count'] != 7599 or result['pending_weight_count'] != 593):
        raise ValueError('Unexpected coverage')
    for pin in [manifest_pin, *audits]:
        load(pin['path'], pin['sha256'])
    if sha256(own.read_bytes()).hexdigest() != own_hash:
        raise RuntimeError('Assembly source changed')
    selection = dict(schema='gamma-partial-cover-v22-explicit-selection-v1',
        construction=result['construction'], target_distance='2/5',
        manifests=[manifest_pin], root_bands=[list(row) for row in ROOT_BANDS],
        selected_audits=audits, selected_event_count=15, selected_weight_count=152,
        assembly_source=dict(path=str(own), sha256=own_hash),
        partial_attempts_not_selected=True, no_weight_interpolation=True,
        numerical_bank_resource_reserve_bytes=750 << 20,
        bounded_final_bookkeeping_reserve_bytes=600 << 20,
        coefficient_replay=False, whole_code_certificate=False)
    sp = write_new(HERE/'hd_field128_gamma30_d40_v22_selection_20261010.json', selection)
    rp = write_new(HERE/'hd_field128_gamma30_d40_partial_cover_v22_20261010.json', result)
    print(json.dumps(dict(selection=sp, index=rp,
        covered=result['covered_weight_count'], pending=result['pending_weight_count'],
        events=result['receipt_count'], union=result['covered_events_union_log2_upper_rational'],
        new_union=result['added_events_union_log2_upper_rational'],
        pending_ranges=result['pending_weight_ranges'])), flush=True)


if __name__ == '__main__':
    main()
