"""Pinned FLOAT gate for the measured bit-syndrome outer / one-state inner.

Reuse only the original unpowered regional curves, never the old outer
bound or the post-refresh construction. This is a complete comparison
for ONE source stratum (r18, arity8, Fmax exactly1), not a whole-code
certificate. The engine and its sources are authenticated before work.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_coupled_state_moment_gate import entry_measure
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power
from hd_two_sum_matrix_moment import ARTIFACT, ARTIFACT_SHA, load_saved


DIRECTORY = Path(__file__).resolve().parent
GROUPS, ARITY, COLUMNS, FLOOR, LIMIT = 18, 8, 16, 1296, 734003
MARK = F(691, 700)
SOURCE_COUNT = 154355220469939768130077702840442839490778956275433472
SOURCE_PINS = {
    'hd_coupled_state_moment_gate.py': 'd2258e48f19bf3c8b7abea6a3d95024c10ce4d7c7286f49d009f52d8a4d7e89a',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
    'hd_two_sum_matrix_moment.py': 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
    'hd_bit_syndrome_kernel.py': 'a27e5557aceaf04273c3e67ce0a01a5d5746c6d01a2ecda8881205ba23541c77',
}
YOUNG_CONTROLS = ((F(1, 2), F(1, 2)), (F(99, 100), F(1, 100)),
                  (F(1, 100), F(99, 100)))


def checked_bank():
    """Validate the reused distribution/scope, not just the array shape."""
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: ' + name)
    data, curves = load_saved()
    required = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), source_count=str(SOURCE_COUNT), source_count_once=True,
        projective_normalizations=1, all_fixed_arity_rowset_vectors=True,
        rowsets_resampled=False, all_regional_state_paths_included=True,
        initial_state=0, terminal_vector=[1, 1], regional_matrix_count=256,
        aggregate_input_weight_floor=FLOOR,
        matrix_entry_order=['00', '01', '10', '11'],
        regional_full_pair_restriction='f<=min(1,floor(m/2))',
        common_potential_used=False, matrix_commutation_assumed=False,
        complete_requested_control=True, directed=False,
        whole_code_certificate=False)
    if any(data.get(key) != value for key, value in required.items()):
        raise ValueError('Saved curves have the wrong one-state distribution or source scope')
    if (curves.shape != (4, GROUPS + 1) or np.isnan(curves).any()
            or np.isposinf(curves).any() or curves[0, 0] != 0):
        raise ValueError('Four complete finite-or-zero occupancy curves required')
    return data, curves


def prefix(log2_x):
    if not np.isfinite(log2_x) or log2_x < 0:
        raise ValueError('Finite activity marker X>=1 required')
    return (log2(SOURCE_COUNT) - LIMIT * log2(float(MARK))
            - (1152 + 1) * log2(255) - FLOOR * log2_x)


def word_bound(matrix):
    power = log_matrix_power(matrix, COLUMNS)
    values = np.logaddexp2.reduce(power[:, 0, :], axis=1)
    return combine_row_bounds(values, GROUPS), values


def run(log2_x, *, engine_sha256, cache_path, information_cover=False,
        max_seconds=120, progress=None):
    """One full DP, plus three cheap predeclared Young-weight recombinations."""
    global_prefix = prefix(log2_x)
    if (type(max_seconds) not in (int, float) or not 0 < max_seconds <= 180
            or type(information_cover) is not bool):
        raise ValueError('Boolean cover and explicit positive budget <=180 seconds required')
    cache_path = Path(cache_path)
    if cache_path.exists():
        raise ValueError('Do not overwrite a frozen branch cache')
    engine_file = DIRECTORY / 'hd_bit_syndrome_gate.py'
    if (len(engine_sha256) != 64 or sha256(engine_file.read_bytes()).hexdigest() != engine_sha256):
        raise ValueError('Reviewed contraction engine changed')
    from hd_bit_syndrome_gate import full_row_bounds, combine_columns
    _, curves = checked_bank()
    started = perf_counter()
    result = dict(construction='bit-syndrome G16 outer, original one-state h4 inner, paired swapped route',
        scope='complete homogeneous r18/arity8/Fmax exactly1 stratum only',
        active_groups=GROUPS, arity=ARITY, columns=COLUMNS, mark=str(MARK),
        log2_activity_marker=float(log2_x), source_count=str(SOURCE_COUNT),
        source_count_once=True, rs_normalizer_exponent=1152, projective_normalizations=1,
        aggregate_input_floor=FLOOR, initial_state=0, terminal_vector=[1, 1],
        rowsets_fixed_across_columns=True, rowset_count=comb(16, 8), extra_rowset_factor=False,
        common_positive_atom_measure=True, branch_sum_before_word_power=True,
        all_boundary_state_paths=True, all_outer_branches=True,
        exact_sparse_restoration=False, information_cover=information_cover,
        raw_regional_banks_built=0, raw_bank_artifact=ARTIFACT,
        raw_bank_artifact_sha256=ARTIFACT_SHA,
        reviewed_source_sha256=dict(SOURCE_PINS, **{'hd_bit_syndrome_gate.py': engine_sha256}),
        prefix_log2=global_prefix, compute_budget_seconds=max_seconds,
        completed_stages=[], complete_requested_control=False,
        directed=False, whole_code_certificate=False)

    def check_time():
        if perf_counter() - started >= max_seconds:
            raise TimeoutError('Fixed bit-syndrome gate exhausted its explicit budget')

    def stage(name):
        result['completed_stages'].append(name)
        if progress:
            progress(dict(stage=name, elapsed_seconds=perf_counter() - started))

    try:
        measure = entry_measure(curves)
        result['state_measure'] = measure
        answer = full_row_bounds(measure, GROUPS, log2_x, lambdas=YOUNG_CONTROLS[0],
            information_cover=information_cover, batch_size=128,
            check_time=check_time, progress=stage)
        if len(answer['rowsets']) != comb(16, 8):
            raise AssertionError('Incomplete labelled row-set coverage')
        branches = answer['branch_column_log2_matrices']
        quiet, nonquiet = branches['quiet'], branches['nonquiet']
        quiet_bound, _ = word_bound(quiet)
        nonquiet_bound, _ = word_bound(nonquiet + ARITY * GROUPS * log2_x)
        controls = []
        best_values = None
        for weights in YOUNG_CONTROLS:
            check_time()
            combined, scales = combine_columns(quiet, nonquiet, GROUPS, ARITY, log2_x, lambdas=weights)
            value, per_rowset = word_bound(combined)
            controls.append(dict(young_weights=list(map(str, weights)),
                young_log2_scales=scales.tolist(), combined_word_log2=value,
                log2_screen=global_prefix + value))
            if best_values is None or value < best_values[0]:
                best_values = value, per_rowset
        with cache_path.open('xb') as handle:
            np.savez_compressed(handle, quiet=quiet, nonquiet=nonquiet,
                rowsets=np.asarray(answer['rowsets']), log2_x=float(log2_x),
                information_cover=information_cover, source_sha256=engine_sha256,
                bank_sha256=ARTIFACT_SHA, output_mark=str(MARK))
        result.update(controls=controls, log2_screen=min(row['log2_screen'] for row in controls),
            quiet_component_log2=global_prefix + quiet_bound,
            nonquiet_component_log2=global_prefix + nonquiet_bound,
            costs=answer['costs'], engine_metadata=answer['metadata'],
            branch_cache=str(cache_path), branch_cache_sha256=sha256(cache_path.read_bytes()).hexdigest(),
            largest_rowsets=[dict(positions=list(answer['rowsets'][int(i)]),
                word_log2=float(best_values[1][i]), root_log2=float(best_values[1][i] / GROUPS))
                for i in np.argsort(best_values[1])[-8:][::-1]],
            complete_requested_control=True)
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter() - started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log2-x', type=float, required=True)
    parser.add_argument('--engine-sha256', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--information-cover', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=120)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Do not overwrite a frozen control or cache')
    answer = run(args.log2_x, engine_sha256=args.engine_sha256, cache_path=args.cache,
        information_cover=args.information_cover, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: answer.get(key) for key in
        ('complete_requested_control', 'log2_screen', 'quiet_component_log2',
         'nonquiet_component_log2', 'elapsed_seconds', 'incomplete_reason')}), flush=True)
