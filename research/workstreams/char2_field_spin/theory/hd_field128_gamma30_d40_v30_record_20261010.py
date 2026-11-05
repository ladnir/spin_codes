"""Publish five pinned selected audits; no remote tree or arithmetic replay."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_selected_index_v30 as selected

HERE = Path(__file__).resolve().parent
HELPER_SHA = '6e914264b8b42c162c93f17d174a979b1b450b99cfa107b4161907e981edc78b'
TEST_SHA = '7d00ada745c8e541170fd2e1266fe79f775c4969dd67141d81386cf16dd12113'
PARENT_SHA = 'ea1900e2f869a258f16404bd34c3882675a78e960e2703a3754820cc4e63737d'
CHECK_SHA = '3a68ec34469cf68f381b4b40cdedf8e1c16220940f14507953845dc23c805056'
RESERVE, OUTPUT_LIMIT, RECORDER_LIMIT = 20 << 30, 512 << 10, 10 << 10
PAIRS = (
    (548, 'b617b50c78ecb735a9aa979a47592a78a194c96dc51ab25b6cb6cd160e926426',
     '867e8602efb5626156664967c62b00ae938fa08d459cd112102968dbe1d59f34'),
    (914, 'aa04947939119bb9ad3ab557b3506412d5ca9131d87135a57bda9c025b790f3d',
     '3e01f37fcb80adad0593d3f5c971768bdceb270919fb2c35b8adc2267869878d'),
    (915, '329aba2e23b6778e12bd54e1f6925b7b517db46031ae9ae6ab456dbd356d6b29',
     '9807e103a1c0ed12222876d53ff046d7eab1598428d2970611c1f9cd31230e85'),
    (916, '61dabb22ec1b587291b96e51cd92f3cd15c18b519893f068fb17500a3d2f527a',
     '15cc7cdfac35f9925e5df5da8790be4e7f725a21a1b50603fbd608de3a90a22b'),
    (1145, 'f94f78e237a5d849adcecdfbf00af72a04a72e0c86c18bb029d57ca56c65a981',
     '84c1efabf5bf6c5a5c7101c9505a13edd1ebd1f0ea16cc88efb6cef72d16a387'),
)


def event_specs():
    return [dict(audit=dict(path=f'remote_receipts_v30/audit_{r}.json', sha256=audit),
                 archive_receipt=dict(path=f'remote_receipts_v30/archive_{r}.json', sha256=archive))
            for r, audit, archive in PAIRS]


def reserve_guard(planned_bytes):
    if shutil.disk_usage(HERE).free-planned_bytes < RESERVE:
        raise OSError('Bookkeeping stopped: preserve 20 GiB after planned writes')


def main():
    reserve_guard(OUTPUT_LIMIT+RECORDER_LIMIT)
    self_path = Path(__file__).resolve()
    self_raw = self_path.read_bytes()
    if len(self_raw) > RECORDER_LIMIT:
        raise OSError('Recorder source exceeds 10 KiB')
    inputs = {self_path:self_raw}
    parent_path = HERE/'hd_field128_gamma30_d40_partial_cover_v29_20261010.json'
    check_path = HERE/'hd_field128_gamma30_d40_partial_cover_v29_20261010_independent_check.json'
    specs = event_specs()
    explicit = [(parent_path,PARENT_SHA), (check_path,CHECK_SHA),
                (Path(selected.__file__).resolve(),HELPER_SHA),
                (HERE/'hd_field128_gamma_selected_index_v30_tests.py',TEST_SHA)]
    explicit += [(HERE/pin['path'],pin['sha256']) for pair in specs for pin in pair.values()]
    for path,digest in explicit:
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned assembly input changed: '+str(path))
        inputs[path] = raw
    parent = json.loads(inputs[parent_path])
    result = selected.extend_v29(parent_path,check_path,specs)
    expected_ranges = [[r,r] for r,unused,unused2 in PAIRS]
    if (result['receipt_count'] != 265 or result['covered_weight_count'] != 8131
            or result['pending_weight_count'] != 61 or result['added_event_count'] != 5
            or result['added_weight_count'] != 5
            or selected.canonical(result['entries'][:260]) != selected.canonical(parent['entries'])
            or [row['weight_range'] for row in result['entries'][260:]] != expected_ranges
            or result['pending_weight_ranges'] != [[521,547],[549,575],[917,919],[1146,1149]]
            or result['covered_events_failure_bits'] != 40
            or result['meets_covered_events_target'] is not True
            or F(result['covered_events_union_log2_upper_rational']) > -40
            or result['whole_code_certificate'] is not False):
        raise ValueError('Unexpected selected events, parent prefix, coverage or global union')
    selection = dict(schema='gamma-v30-explicit-selected-remote-evidence-selection-v1',
        parent=dict(path=parent_path.name,sha256=PARENT_SHA),
        parent_independent_check=dict(path=check_path.name,sha256=CHECK_SHA),
        validator=dict(path=Path(selected.__file__).name,sha256=HELPER_SHA),
        validator_tests=dict(path='hd_field128_gamma_selected_index_v30_tests.py',sha256=TEST_SHA),
        recorder=dict(path=self_path.name,sha256=sha256(self_raw).hexdigest()),
        event_specs=specs,added_events=5,added_weights=5,
        previous_260_records_and_order_unchanged=True,selected_weights=[r for r,unused,unused2 in PAIRS],
        actual_ceiling_union_failure_bits=40,individual_requested_targets_preserved=True,
        composite_bands_included=False,nonoverlapping_weight_ranges=True,
        mathematical_audits_inherited_from_remote=True,archive_verification_inherited_from_remote=True,
        historical_origin_verified=False,local_tree_replay=False,local_archive_revalidation=False,
        coefficient_replay=False,normalized_regional_replay=False,transfer_recomputed=False,
        whole_code_certificate=False,combined_selection_index_byte_limit=OUTPUT_LIMIT,
        recorder_source_byte_limit=RECORDER_LIMIT,minimum_free_bytes_after_planned_writes=RESERVE,
        storage_scope='Compact records only; remote audits and durable archive verification are inherited')
    outputs = [(HERE/'hd_field128_gamma30_d40_v30_selection_20261010.json',selection),
               (HERE/'hd_field128_gamma30_d40_partial_cover_v30_20261010.json',result)]
    encoded = [(path,(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()) for path,value in outputs]
    total = sum(len(raw) for unused,raw in encoded)
    if total > OUTPUT_LIMIT:
        raise OSError('Selection and index exceed 512 KiB')
    reserve_guard(total)
    if any(path.exists() for path,unused in encoded):
        raise FileExistsError('Refusing to replace immutable selection or index')
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input changed during assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    remaining = total
    for path,raw in encoded:
        reserve_guard(remaining)
        with path.open('xb') as stream:
            stream.write(raw)
        remaining -= len(raw)
        if path.read_bytes() != raw:
            raise RuntimeError('Published bytes differ: '+str(path))
        print(json.dumps(dict(path=str(path),bytes=len(raw),sha256=sha256(raw).hexdigest())),flush=True)
    if any(path.read_bytes() != raw for path,raw in inputs.items()):
        raise RuntimeError('Input changed before closing assembly')
    selected.frozen.producer.check_source_pins(result['aggregation_source_sha256'])
    print(json.dumps(dict(selection_index_bytes=total,recorder_bytes=len(self_raw),
        **{key:result[key] for key in ('receipt_count','covered_weight_count','pending_weight_count',
            'pending_weight_ranges','covered_events_union_log2_upper_rational',
            'added_events_union_log2_upper_rational','whole_code_certificate')})),flush=True)


if __name__ == '__main__':
    main()
