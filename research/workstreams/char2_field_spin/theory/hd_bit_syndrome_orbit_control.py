"""One pinned per-orbit FLOAT control, with the original one-state provider.

This is a selected r18/A8/Fmax1 class. No native construction is changed,
and no whole-code or directed certificate is claimed. The permitted marks
are the saved control, one alternative, and one subsequent interior mark
motivated by the alternative's dominant all-live path; no optimizer.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy

from hd_bit_syndrome_gate_run import (DIRECTORY, GROUPS, ARITY, COLUMNS, FLOOR,
    LIMIT, SOURCE_COUNT, YOUNG_CONTROLS, checked_bank, word_bound)
from hd_coupled_state_moment_gate import entry_measure, MomentFitError
from hd_two_sum_gf16_factorized import forced_nonzero_regions
from hd_two_sum_late_alpha_source_f1 import restrict_regions, load_source_count


SOURCE_PINS = {
    'hd_bit_syndrome_gate_run.py': '4b938c7d8fbd46c81ffb4fc821143c88b95722e12d1228b442a23d3583b00a23',
    'hd_bit_syndrome_gate.py': '3f8fa49998b1aacc00a964e5c6b83611d62510424afeca98875979aea8063505',
    'hd_two_sum_gf16_factorized.py': 'd34f054a6a1a07cec5ae97c75acdf475763a5c5e5ccab85d8be05cdd281e5da2',
    'hd_two_sum_late_alpha_source_f1.py': '77d3f7da4f7439766d0bb18b53cebcd0960697094c0e30361b064884b5ee1987',
    'hd_xor_multigroup.py': 'e64a98c2db2a04d34ee61a7814825a85fba6724d281e0963215c6e59442d1bd0',
    'hd_shifted_swap_sparse.py': '9b8875ec1433a2bc79dbd1e42d3b90815528a70ce0234069cfd03d61b030c58c',
    'hd_bucket_feedback.py': '66a1e97608046ae9e1fa6e37f6b6e5634de586fd03ae819278014d6a8ae16d38',
    'hd_common_refresh_screen.py': 'e1af79f18c7634dd3c2489185bb6c95d0b54a207d3a360fcfeff17b94093922f',
    'hd_common_refresh_principal.py': '1ecd8f4243b81ca095f045d6c6b9c28529900bcfff1c5fa4382a2cd8ca259cb5',
    'hd_grouped_quiet.py': '4d84222bcaaeab1d134bb3ab8235d4082f39ebb9970e47efa110f52f9bf5a369',
    'hd_grouped_subfield_screen.py': '74106dec2416e5f41bf07a80d922796e6cba659aa8e491ee0cc0452f729b5964',
    'hd_xor_column_probe.py': '06b3ecbd39c2aa917851745fa4f197f5c5641ea266d6700362a39fb8ead1b563',
    'hd_systematic_mds.py': 'fe04df5613156b79215784e960e8a85e3cf6d4e17f8dc0c268c0eb911e66c5e2',
    'hd_blockwise_shared_transfer.py': '0245af670491b5aaab89cb6d9e620a2ebc490907eef4dea8dc1ee804085a260b',
    'hd_blockwise_maps.py': '1f80f1f8c6f89509c69043b4542ba9598d56d6e1afe2805159c2c50848becab4',
    'hd_sumrank_maps.py': '704713311f17fe327aeeece946326ed53e2487e752f86fe96aeb82e135d2577b',
}


def global_prefix(mark, log2_x):
    mark = F(mark)
    if not 0 < mark < 1 or not np.isfinite(log2_x) or log2_x < 0:
        raise ValueError('A rational mark in (0,1) and finite X>=1 required')
    return (log2(SOURCE_COUNT) - LIMIT * log2(float(mark))
            - 1153 * log2(255) - FLOOR * log2_x)


def run(mark, log2_x, *, engine_sha256, cache_path, max_seconds=180, progress=None):
    mark = F(mark)
    prefix = global_prefix(mark, log2_x)
    if mark not in (F(691, 700), F(497, 500), F(993, 1000)):
        raise ValueError('Only a saved, prescribed alternative, or interior diagnostic mark is permitted')
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 180:
        raise ValueError('Explicit positive budget <=180 seconds required')
    cache_path = Path(cache_path)
    if cache_path.exists():
        raise ValueError('Do not overwrite a frozen cache')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: ' + name)
    engine = DIRECTORY / 'hd_bit_syndrome_orbit_gate.py'
    if len(engine_sha256) != 64 or sha256(engine.read_bytes()).hexdigest() != engine_sha256:
        raise ValueError('Reviewed per-orbit engine changed')
    from hd_bit_syndrome_orbit_gate import full_row_bounds
    from hd_bit_syndrome_gate import combine_columns
    # Also verifies the frozen original distribution and underlying helper pins.
    _, saved_curves = checked_bank()
    if load_source_count() != SOURCE_COUNT:
        raise ValueError('Exact source partition changed')
    started = perf_counter()
    result = dict(construction='bit-syndrome G16 outer, original one-state h4 inner, paired swapped route',
        scope='complete homogeneous r18/arity8/Fmax exactly1 stratum only',
        active_groups=GROUPS, arity=ARITY, columns=COLUMNS, mark=str(mark),
        log2_activity_marker=float(log2_x), field_order=256, forbidden_weight=LIMIT,
        message_length=1048576, codeword_length=2097152, record_bytes=16,
        inner='h4: packet64, state16, four RS[16,4] maps, independent D2 H16 D1 update',
        route='alternating even/odd adjacent pairs; independent pair permutations and swap bits',
        regional_state_basis=['zero', 'normalized uniform nonzero'],
        regional_matrix_count=256, packets_per_region=128,
        provider_binary64_mark_exact=str(F(float(mark))), zero_entries_are_none=True,
        runtime_versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        source_count=str(SOURCE_COUNT), source_count_once=True, rs_normalizer_exponent=1152,
        projective_normalizations=1, aggregate_input_floor=FLOOR,
        initial_state=0, terminal_vector=[1, 1], rowset_count=comb(16, 8),
        rowsets_fixed_across_columns=True, chronological_rowsets_collapsed=False,
        extra_rowset_factor=False, common_positive_atom_measure=True,
        branch_sum_before_word_power=True, all_outer_branches=True,
        all_boundary_state_paths=True, exact_sparse_restoration=False,
        information_cover=True, per_orbit_coefficient_caps=True,
        regional_alpha=1, regional_full_pair_restriction='f<=min(1,floor(m/2))',
        matrix_entry_order=['00', '01', '10', '11'], common_potential_used=False,
        matrix_commutation_assumed=False, prefix_log2=prefix,
        reviewed_source_sha256=dict(SOURCE_PINS, **{'hd_bit_syndrome_orbit_gate.py': engine_sha256}),
        compute_budget_seconds=max_seconds, completed_stages=[],
        raw_regional_banks_built=0, complete_requested_control=False,
        directed=False, whole_code_certificate=False)

    def check_time():
        if perf_counter() - started >= max_seconds:
            raise TimeoutError('Per-orbit control exhausted its explicit budget')

    def stage(name):
        result['completed_stages'].append(name)
        if progress:
            progress(dict(stage=name, elapsed_seconds=perf_counter() - started))

    try:
        curves = saved_curves
        if mark != F(691, 700):
            regions = restrict_regions(forced_nonzero_regions(mark, GROUPS))
            curves = np.array([regions[m, 0].ravel() for m in range(GROUPS + 1)]).T
            result['raw_regional_banks_built'] = 1
        if curves.shape != (4, 19) or np.isnan(curves).any() or np.isposinf(curves).any():
            raise ValueError('Complete finite-or-zero four-entry bank required')
        result['regional_log2_curves'] = [[None if x == -np.inf else float(x) for x in row] for row in curves]
        stage('original_one_state_bank_complete')
        check_time()
        measure = entry_measure(curves)
        result['state_measure'] = measure
        answer = full_row_bounds(measure, GROUPS, log2_x, information_cover=True,
            lambdas=YOUNG_CONTROLS[0], batch_size=128, check_time=check_time, progress=stage)
        if len(answer['rowsets']) != comb(16, 8):
            raise AssertionError('Incomplete labelled source-set coverage')
        quiet = answer['branch_column_log2_matrices']['quiet']
        nonquiet = answer['branch_column_log2_matrices']['nonquiet']
        quiet_bound, _ = word_bound(quiet)
        nonquiet_bound, _ = word_bound(nonquiet + ARITY * GROUPS * log2_x)
        controls, best_values = [], None
        for weights in YOUNG_CONTROLS:
            check_time()
            combined, scales = combine_columns(quiet, nonquiet, GROUPS, ARITY, log2_x, lambdas=weights)
            value, values = word_bound(combined)
            controls.append(dict(young_weights=list(map(str, weights)),
                young_log2_scales=scales.tolist(), combined_word_log2=value, log2_screen=prefix + value))
            if best_values is None or value < best_values[0]:
                best_values = value, values
        with cache_path.open('xb') as handle:
            np.savez_compressed(handle, quiet=quiet, nonquiet=nonquiet, rowsets=np.asarray(answer['rowsets']),
                log2_x=float(log2_x), information_cover=True, per_orbit_caps=True,
                engine_sha256=engine_sha256, output_mark=str(mark), regional_log2_curves=curves)
        result.update(controls=controls, log2_screen=min(row['log2_screen'] for row in controls),
            quiet_component_log2=prefix + quiet_bound, nonquiet_component_log2=prefix + nonquiet_bound,
            costs=answer['costs'], engine_metadata=answer['metadata'],
            branch_cache=str(cache_path), branch_cache_sha256=sha256(cache_path.read_bytes()).hexdigest(),
            largest_rowsets=[dict(positions=list(answer['rowsets'][int(i)]),
                word_log2=float(best_values[1][i]), root_log2=float(best_values[1][i] / GROUPS))
                for i in np.argsort(best_values[1])[-8:][::-1]], complete_requested_control=True)
    except (TimeoutError, MomentFitError) as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter() - started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mark', type=F, required=True)
    parser.add_argument('--log2-x', type=float, required=True)
    parser.add_argument('--engine-sha256', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Do not overwrite a frozen control or cache')
    answer = run(args.mark, args.log2_x, engine_sha256=args.engine_sha256, cache_path=args.cache,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: answer.get(key) for key in ('complete_requested_control', 'log2_screen',
        'quiet_component_log2', 'nonquiet_component_log2', 'elapsed_seconds', 'incomplete_reason')}), flush=True)
