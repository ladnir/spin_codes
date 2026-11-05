"""Saved eight-cone bank -> fixed live compression -> ordered moment gate.

Only one predeclared live reference and the prescribed 00(5,18), 01(2,18),
10(2,18), constant11 moment laws are used. The frozen positive polynomial
engine retains all masks and boundary paths, Q(X)^16 and degree>=1296.
This FLOAT screen builds no packet/regional bank and transfers no distance
receipt. A raw-bank hash must be supplied explicitly by the caller.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_post_refresh_compress import (STATE_ORDER, FIXED_LOG_WEIGHTS, checked_bank,
    compress_bank, entrywise_curves)
from hd_coupled_state_moment_gate import MomentFitError, entry_measure, polynomial_control
from hd_two_sum_gf16_mean_column import raw_holder_weights
from hd_two_sum_matrix_moment import load_saved as load_source_count
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
CONSTRUCTION = 'proposed-g16-post-refresh-coupled-inner'
GROUPS, ARITY, COLUMNS, FLOOR = 18, 8, 16, 1296
MARK = Fraction(691, 700)
SOURCE_COUNT = 154355220469939768130077702840442839490778956275433472
PROVIDER_SHA = '5942384325e38e77cf50db1a02b4b8ca30e50002256d9ceac25359158b2615ee'
PROVIDER_DEPENDENCIES = {
    'hd_post_refresh_packet.py': 'b71bc026e32fb673478958df2fc0ec243c704024f796988d773d30c4fb607fbd',
    'hd_coupled_state_packet.py': 'faab3a80e582f0c1864b79f63b256788d99a5034f15bfbd9e781940b43bc2a43',
    'hd_coupled_state_regions.py': 'd843965b9e74d1aa86b80d09c7ab4e3a0afa594dc2c91057a7fa7497c97fa01f',
}
SOURCE_PINS = dict(PROVIDER_DEPENDENCIES, **{
    'hd_post_refresh_regions.py': PROVIDER_SHA,
    'hd_post_refresh_compress.py': '411271e8a507021ae81e030d31a1f20ea611da3c96d976543b8c389e1fb20fcc',
    'hd_coupled_state_moment_gate.py': 'd2258e48f19bf3c8b7abea6a3d95024c10ce4d7c7286f49d009f52d8a4d7e89a',
    'hd_two_sum_gf16_mean_column.py': 'fa6377002fa52fbfcb980ddd4b72d30df16b6c12c7c6db24eed54554e61a1593',
    'hd_two_sum_matrix_moment.py': 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67',
    'hd_two_sum_gf16_state_envelope.py': '5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7',
    'hd_two_sum_late_alpha_quiet.py': '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
})


def decode_saved_bank(data):
    """Fail closed on construction, complete coefficient domain, or provenance."""
    required = dict(construction=CONSTRUCTION, mark=str(MARK), maximum_regional_occupancy=18,
        full_pair_maximum=1, regional_alpha=1, state_order=list(STATE_ORDER), initial_state='Z',
        terminal_vector=[1]*8, regional_histogram_count=36, zero_matrix_entries_encoded_as_none=True,
        raw_packet_banks_built=1, parameter_fits=0, complete_requested_bank=True,
        directed=False, whole_code_certificate=False, source_sha256=PROVIDER_SHA,
        reviewed_source_sha256=PROVIDER_DEPENDENCIES)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Exact complete pinned post-refresh r18/f<=1 raw bank required')
    metadata = data.get('provider_metadata', {})
    expected_metadata = dict(construction=CONSTRUCTION, state_order=list(STATE_ORDER),
        regional_alpha=1, total_packets=128, total_pairs=4096, global_normalization_count=1,
        coefficient_maximum=18, coefficient_full_max=1, orientation_averaged_in_packet=True,
        source_geometry_resampled=False, states_persist_across_regions=True,
        cone_P_arbitrary_diagonal_probability=True, cone_G_arbitrary_first_nonzero_probability=True,
        cone_F_conditional_second_flat_density=True, probability_rows_assumed=False,
        directed=False, whole_code_certificate=False)
    if any(metadata.get(key) != value for key, value in expected_metadata.items()):
        raise ValueError('Raw provider chronology or cone metadata changed')
    bank = {}
    for record in data['regional_log2_matrices']:
        key = record['active_weight'], record['full_pairs']
        if any(type(value) is not int for value in key) or key in bank:
            raise ValueError('Unique integer histogram keys required')
        bank[key] = np.array([[-np.inf if value is None else value for value in row]
                              for row in record['log2_matrix']], dtype=float)
    return checked_bank(bank)


def load_saved_bank(path, expected_sha256):
    if (not isinstance(expected_sha256, str) or len(expected_sha256) != 64
            or any(c not in '0123456789abcdef' for c in expected_sha256)):
        raise ValueError('Explicit lowercase SHA256 of the frozen raw bank required')
    path = Path(path)
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Saved raw bank does not match its requested pin')
    data = json.loads(raw)
    return data, decode_saved_bank(data)


def check_source_pins():
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)


def run(bank_path, bank_sha256, *, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive resource budget at most thirty seconds required')
    check_source_pins()
    raw_data, bank = load_saved_bank(bank_path, bank_sha256)
    count_data, _ = load_source_count()
    if int(count_data['source_count']) != SOURCE_COUNT:
        raise ValueError('Pinned exact Fmax=1 source count changed')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Saved post-refresh moment gate exhausted its budget')
    prefix = log2(SOURCE_COUNT)-LIMIT*log2(float(MARK))-log2(255)-1152*log2(255)
    result = dict(construction=CONSTRUCTION, scope='complete homogeneous r18/arity8/Fmax1 class only',
        active_groups=GROUPS, arity=ARITY, source_full_pair_max=1, source_count=str(SOURCE_COUNT),
        source_count_condition='Fmax exactly1', mark=str(MARK), alpha=1, regional_alpha=1,
        state_order=list(STATE_ORDER), fixed_live_log_weights=list(FIXED_LOG_WEIGHTS),
        weights_selected_before_raw_bank=True, common_live_reference_all_regions=True,
        moment_secants={'00': [5, 18], '01': [2, 18], '10': [2, 18], '11': 'constant'},
        fixed_rowsets_across_columns=True, raw_holder_exponent=16, extra_rowset_factor=False,
        source_count_once=True, projective_normalizations=1, rs_normalizer_exponent=1152,
        aggregate_surrogate_floor=FLOOR, all_outer_masks=True, all_boundary_state_paths=True,
        terminal_vector=[1, 1], extra_terminal_factor_log2=0., prefix_log2=prefix,
        offdiagonal_empty_zeros_retained=False, parameter_fits=0, raw_regional_banks_built=0,
        old_construction_receipts_transferred=False, reviewed_source_sha256=SOURCE_PINS,
        raw_bank_artifact=str(Path(bank_path)), raw_bank_artifact_sha256=bank_sha256,
        raw_bank_source_sha256=raw_data['source_sha256'], regional_histogram_count=len(bank),
        compute_budget_seconds=max_seconds, complete_requested_control=False,
        directed=False, whole_code_certificate=False)
    try:
        check_time()
        compressed = compress_bank(bank)
        curves = entrywise_curves(compressed['regions'])
        result.update(compression=compressed['metadata'],
            compressed_entry_log2_curves=[[None if value == -np.inf else float(value)
                                          for value in row] for row in curves])
        measure = entry_measure(curves)
        result['state_measure'] = measure
        check_time()
        control = polynomial_control(raw_holder_weights(), measure, prefix, check_time=check_time)
        result.update(control)
        check_time()
        result['complete_requested_control'] = True
        if progress:
            progress(dict(stage='post_refresh_gate_complete', log2_screen=result['log2_screen']))
    except (TimeoutError, MomentFitError) as error:
        result['incomplete_reason'] = str(error)
        result['incomplete_kind'] = type(error).__name__
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank', type=Path, required=True)
    parser.add_argument('--bank-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen gate')
    answer = run(args.bank, args.bank_sha256, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
