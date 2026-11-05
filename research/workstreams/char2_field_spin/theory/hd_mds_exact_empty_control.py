"""One fresh same-mark mixed7/8/9 control with the exact empty-step cap.

No encoder change, cone change, witness search, or imported distance bound.
All kernels use one freshly checked common matrix-valued atom measure.
This is a selected-class FLOAT control, not a whole-code certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_activity_control import constant_columns, marked_columns, norm
from hd_bit_syndrome_arity7_gate import full_row_bounds as arity7
from hd_bit_syndrome_arity9_gate import full_row_bounds as arity9
from hd_bit_syndrome_directed_recheck import independent_source_counts
from hd_bit_syndrome_mixed_control import PINS
from hd_bit_syndrome_orbit_gate import full_row_bounds as arity8
from hd_coupled_state_moment_gate import entry_measure
from hd_mds_exact_empty_regions import forced_nonzero_regions
from hd_two_sum_late_alpha_source_f1 import restrict_regions

DIRECTORY = Path(__file__).resolve().parent
MARK, GROUPS, LIMIT = F(9927, 10000), 18, 734003
EXTRA_PINS = {
    'hd_bit_syndrome_arity7_gate.py': 'b20e009797cc6a0986dbd594ff67a061197c1867c08eecf2922aa8085d1a085f',
    'hd_bit_syndrome_counted_general.py': '980ec6e2d20f7eae25a2819396cbd2b2a5f901560cba4cf5dfeee20267de6283',
    'hd_bit_syndrome_activity_control.py': '533fe489bacad9af0505dc1d17bef94e4a0d0a7807170ad33a51e036417346b3',
    'hd_mds_empty_weight.py': '6c76d7513942d67d10ebe7b5610517525fdb29f5de42e416db474b92b49e3920',
    'hd_xor_state_local.py': '61e90b22736b50adb2ad7ff4e8e61fd09df592297b39a2b5220eb160abd0b502',
}


def run(cache_path, *, adapter_sha, max_seconds=120, progress=None):
    cache_path = Path(cache_path)
    if cache_path.exists() or not 0 < max_seconds <= 180:
        raise ValueError('New cache and explicit bounded budget required')
    pins = dict(PINS, **EXTRA_PINS, **{'hd_mds_exact_empty_regions.py': adapter_sha})
    for name, expected in pins.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)
    started, stages = perf_counter(), []
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Exact-empty selected control exceeded its budget')
    def stage(label):
        row = dict(stage=label, elapsed_seconds=perf_counter()-started)
        stages.append(row)
        if progress:
            progress(row)

    regions, provider = forced_nonzero_regions(MARK, GROUPS, check_time=check_time)
    regions = restrict_regions(regions)
    curves = np.array([regions[m, 0].ravel() for m in range(GROUPS+1)]).T
    measure = entry_measure(curves)
    stage('fresh_common_measure_exact_empty_packet')
    seven = arity7(measure, GROUPS, check_time=check_time)
    stage('all_labelled_arity7')
    eight = arity8(measure, GROUPS, 0., information_cover=True,
        lambdas=(F(1, 8), F(7, 8)), check_time=check_time)
    stage('all_labelled_arity8')
    nine = arity9(measure, GROUPS, check_time=check_time)
    stage('all_labelled_arity9')
    columns = {7: seven['column_log2_matrices'], 8: eight['column_log2_matrices'],
               9: nine['column_log2_matrices']}
    source_count = independent_source_counts()['F_le1']
    prefix = log2(source_count)-LIMIT*log2(MARK)-log2(255)
    base, _ = constant_columns(measure)
    controls = []
    for marker in (F(1), F(5, 2)):
        values = {a: norm(marked_columns(column, base, marker), a, marker)
                  if a == 7 else norm(column, a) for a, column in columns.items()}
        controls.append(dict(arity7_marker=str(marker), arity8_marker='1', arity9_marker='1',
            per_arity_log2={str(a): prefix+v for a, v in values.items()},
            mixed_log2=prefix+GROUPS*np.logaddexp2.reduce([v/GROUPS for v in values.values()])))
    with cache_path.open('xb') as handle:
        np.savez_compressed(handle, column7=columns[7], column8=columns[8], column9=columns[9],
            quiet8=eight['branch_column_log2_matrices']['quiet'],
            nonquiet8=eight['branch_column_log2_matrices']['nonquiet'],
            rowsets7=np.asarray(seven['rowsets']), rowsets8=np.asarray(eight['rowsets']),
            rowsets9=np.asarray(nine['rowsets']))
    return dict(construction=seven['metadata']['outer_construction'], groups=GROUPS,
        arities=[7, 8, 9], source_class='Fmax<=1', source_count=str(source_count),
        mark=str(MARK), input_marker=1, state_measure=measure,
        regional_log2_curves=[[None if x == -np.inf else float(x) for x in row] for row in curves],
        common_measure=True, common_prefix_log2=prefix, controls=controls,
        empty_packet_provider=provider, stages=stages,
        elapsed_seconds=perf_counter()-started, reviewed_sources=pins,
        cache=str(cache_path), cache_sha256=sha256(cache_path.read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        output_mark_fits=0, native_encoder_changed=False, directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--adapter-sha256', required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Do not overwrite a saved control')
    result = run(args.cache, adapter_sha=args.adapter_sha256,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('controls', 'elapsed_seconds')}), flush=True)
