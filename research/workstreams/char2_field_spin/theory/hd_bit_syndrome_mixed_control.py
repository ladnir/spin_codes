"""One interior output-marker control between the a8/a9 obstructions.

The124/125 predecessor is preserved in the adjacent frozen snapshot.
It cleared a9 but made the a8 nonquiet branch fail. Fresh original-one-state
curves at9927/10000 now supply ONE common FLOAT atom
measure to both arities. Two prescribed a8 Young choices reuse its branch
matrices. Arity9 keeps the full actual-fiber bound, including outside zero.
This is exploratory numerical evidence, never a certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_arity9_gate import full_row_bounds as arity9
from hd_bit_syndrome_directed_recheck import independent_source_counts
from hd_bit_syndrome_gate import combine_columns
from hd_bit_syndrome_orbit_control import SOURCE_PINS
from hd_bit_syndrome_orbit_gate import full_row_bounds as arity8
from hd_coupled_state_moment_gate import entry_measure
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power
from hd_two_sum_gf16_factorized import forced_nonzero_regions
from hd_two_sum_late_alpha_source_f1 import restrict_regions

DIRECTORY = Path(__file__).resolve().parent
MARK, GROUPS, LIMIT = F(9927, 10000), 18, 734003
WEIGHTS = ((F(1, 8), F(7, 8)), (F(1, 10), F(9, 10)))
PINS = dict(SOURCE_PINS, **{
    'hd_coupled_state_moment_gate.py': 'd2258e48f19bf3c8b7abea6a3d95024c10ce4d7c7286f49d009f52d8a4d7e89a',
    'hd_bit_syndrome_directed_recheck.py': '224d1bb24dd5a90f4d141c68349472265f66ec1b2270a1d8a5d3c93a3780eadd',
    'hd_bit_syndrome_arity9_gate.py': 'a02369609ee860cec6fd236e952e82e841ef776cb369c4735a88fbb93d186e48',
    'hd_bit_syndrome_cross_fibers.py': 'd8658b8c3a67d2925e776b63bd62b637456ea9e22c8c75d64dda1ba958e9104e',
    'hd_bit_syndrome_orbit_gate.py': '1aad0e0c7b88397f18f5c81806ba04b4a2aab64015396e740af4155ea7636023',
    'hd_bit_syndrome_kernel.py': 'a27e5557aceaf04273c3e67ce0a01a5d5746c6d01a2ecda8881205ba23541c77',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
})


def run(cache_path, *, max_seconds=180, progress=None):
    cache_path = Path(cache_path)
    if cache_path.exists() or not 0 < max_seconds <= 300:
        raise ValueError('A new cache and bounded run budget required')
    for name, expected in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)
    started, stages = perf_counter(), []

    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Mixed-arity control exhausted its explicit budget')

    def stage(name):
        row = dict(stage=name, elapsed_seconds=perf_counter()-started)
        stages.append(row)
        if progress:
            progress(row)

    regions = restrict_regions(forced_nonzero_regions(MARK, GROUPS))
    curves = np.array([regions[m, 0].ravel() for m in range(GROUPS+1)]).T
    if curves.shape != (4, 19) or np.isnan(curves).any() or np.isposinf(curves).any():
        raise ValueError('Complete original nonzero-amplitude curves required')
    measure = entry_measure(curves)
    stage('fresh_common_original_inner_measure')
    a8 = arity8(measure, GROUPS, 0., information_cover=True, lambdas=WEIGHTS[0],
                check_time=check_time, progress=stage)
    a9 = arity9(measure, GROUPS, check_time=check_time)
    stage('all_labelled_arity9_full_fibers')
    counts = independent_source_counts()
    prefix = log2(counts['F_le1'])-LIMIT*log2(float(MARK))-log2(255)
    norm9 = a9['combined_word_log2']-72*GROUPS*log2(255)
    controls = []
    for weights in WEIGHTS:
        combined, _ = combine_columns(a8['branch_column_log2_matrices']['quiet'],
            a8['branch_column_log2_matrices']['nonquiet'], GROUPS, 8, 0., lambdas=weights)
        word = log_matrix_power(combined, 16)
        norm8 = combine_row_bounds(np.logaddexp2.reduce(word[:, 0, :], axis=1), GROUPS)-64*GROUPS*log2(255)
        mixed = GROUPS*np.logaddexp2(norm8/GROUPS, norm9/GROUPS)
        controls.append(dict(arity8_young=list(map(str, weights)), arity8_log2=prefix+norm8,
                             mixed_arity8_or9_log2=prefix+mixed))
    with cache_path.open('xb') as handle:
        np.savez_compressed(handle, quiet8=a8['branch_column_log2_matrices']['quiet'],
            nonquiet8=a8['branch_column_log2_matrices']['nonquiet'],
            column9=a9['column_log2_matrices'], words9=a9['per_rowset_word_log2'],
            rowsets8=np.asarray(a8['rowsets']), rowsets9=np.asarray(a9['rowsets']))
    return dict(construction=a9['metadata']['outer_construction'], groups=GROUPS,
        source_class='Fmax<=1', mark=str(MARK), input_marker=1, arity9_log2=prefix+norm9,
        controls=controls, state_measure=measure,
        regional_log2_curves=[[None if x == -np.inf else float(x) for x in row] for row in curves],
        common_measure_for_both_arities=True, source_count_once=True, projective_normalizations=1,
        mixed_method='18th power of sum of normalized per-arity18th roots',
        stages=stages, elapsed_seconds=perf_counter()-started,
        a8_metadata=a8['metadata'], a9_metadata=a9['metadata'], reviewed_sources=PINS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), cache=str(cache_path),
        cache_sha256=sha256(cache_path.read_bytes()).hexdigest(), directed=False,
        whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.cache.exists() or args.output.exists():
        parser.error('Do not overwrite a saved control')
    result = run(args.cache, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('arity9_log2', 'controls', 'elapsed_seconds')}), flush=True)
