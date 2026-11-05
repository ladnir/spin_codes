"""Assemble exactly the 36 audited additions selected for partial index v20."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_selected_index_v20 as index


HERE = Path(__file__).resolve().parent
MANIFESTS = (
    ('hd_field128_gamma30_d40_kernel_v20_band_manifest_20261010.json',
     'caafde6960733825d132e1694c1a6a14d9891a7e95c76c84c6f4bb9df43e558b', 'events', 15),
    ('hd_field128_gamma30_d40_upper_gap_5984_6677_manifest_20261010.json',
     'ed45de28596d630ff0b6937e34d16a14eb61667bae50555a58dd4911e06b0d26',
     'individually_audited_events', 11),
)
ROOT_BANDS = ((3424, 3487), (3488, 3551), (3552, 3615), (3616, 3647),
              (3648, 3711), (3712, 3775), (3776, 3839), (3840, 3903),
              (5280, 5343), (5344, 5344))


def load(path, expected=None):
    path = Path(path).resolve()
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if expected is not None and digest != expected:
        raise RuntimeError('Immutable input changed: '+str(path))
    return json.loads(raw), dict(path=str(path), sha256=digest)


def main():
    manifests, audits = [], []
    for name, digest, key, count in MANIFESTS:
        manifest, pin = load(HERE/name, digest)
        if len(manifest[key]) != count:
            raise ValueError('Unexpected selected event count')
        manifests.append(pin)
        for event in manifest[key]:
            unused, audit_pin = load(event['independent_audit']['path'],
                                     event['independent_audit']['sha256'])
            audits.append(audit_pin)
    for lo, hi in ROOT_BANDS:
        folder = HERE/f'hd_field128_gamma30_d40_weights{lo}_{hi}_midpoint_band_20261010'
        audit, pin = load(folder/'independent_audit.json')
        if audit['requested_weight_range'] != [lo, hi]:
            raise ValueError('Unexpected root audit interval')
        audits.append(pin)
    if len(audits) != 36 or len({p['path'] for p in audits}) != 36:
        raise ValueError('Selection is not 36 distinct audited events')
    result = index.extend_v19(HERE/'hd_field128_gamma30_d40_partial_cover_v19_20261010.json',
                             [p['path'] for p in audits])
    if result['added_event_count'] != 36 or result['added_weight_count'] != 1731:
        raise ValueError('Unexpected new coverage')
    for pin in [*manifests, *audits]:
        load(pin['path'], pin['sha256'])
    selection = dict(schema='gamma-partial-cover-v20-explicit-selection-v1',
        construction=result['construction'], target_distance='2/5',
        manifests=manifests, root_bands=[list(pair) for pair in ROOT_BANDS],
        selected_audits=audits, selected_event_count=36, selected_weight_count=1731,
        assembly_source=dict(path=str(Path(__file__).resolve()),
                             sha256=sha256(Path(__file__).read_bytes()).hexdigest()),
        excluded_root_partial_interval=[3904, 3936],
        no_weight_interpolation=True, coefficient_replay=False, whole_code_certificate=False)
    selection_pin = index.frozen.producer.profiles.write_new(
        HERE/'hd_field128_gamma30_d40_v20_selection_20261010.json', selection)
    result_pin = index.frozen.producer.profiles.write_new(
        HERE/'hd_field128_gamma30_d40_partial_cover_v20_20261010.json', result)
    print(json.dumps(dict(selection=selection_pin, index=result_pin,
        covered=result['covered_weight_count'], pending=result['pending_weight_count'],
        events=result['receipt_count'], union=result['covered_events_union_log2_upper_rational'],
        pending_ranges=result['pending_weight_ranges'])), flush=True)


if __name__ == '__main__':
    main()
