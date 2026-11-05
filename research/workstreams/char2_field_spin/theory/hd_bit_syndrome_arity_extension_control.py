"""Bounded extra-arity FLOAT control under the frozen mixed8/9 measure.

Every comparison uses r18,F<=1,z9927/10000,X1. No provider or refit runs.
The mixed sum is normalized per group before the eighteenth power, not a
sum of homogeneous bounds. All results remain partial FLOAT evidence.
"""
from fractions import Fraction as F
from hashlib import sha256
import importlib
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_directed_recheck import independent_source_counts
from hd_bit_syndrome_gate import combine_columns
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power

DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_bit_syndrome_mixed89_9927_x1_float.json'
PROPOSAL_SHA = 'd94ccf2b9a9bc073b67c90256608fa67ef7be02f44efc47eede301f627c4e2c3'
CACHE_SHA = '0cdf0b60933ae1696c434afb4986caa9bf29f986a2f520762ba5245b1f0f6067'
COUNTED_SHA = '980ec6e2d20f7eae25a2819396cbd2b2a5f901560cba4cf5dfeee20267de6283'
GROUPS, LIMIT, MARK = 18, 734003, F(9927, 10000)


def run(arity, engine_sha, cache_path, *, max_seconds=180, progress=None):
    if arity != 7 and not 10 <= arity <= 16:
        raise ValueError('This adapter extends a7 or a10..16 only')
    name = 'hd_bit_syndrome_arity7_gate' if arity == 7 else 'hd_bit_syndrome_higharity_gate'
    if len(engine_sha) != 64 or sha256((DIRECTORY/(name+'.py')).read_bytes()).hexdigest() != engine_sha:
        raise ValueError('Reviewed extra-arity engine changed')
    if sha256((DIRECTORY/PROPOSAL_NAME).read_bytes()).hexdigest() != PROPOSAL_SHA:
        raise ValueError('Frozen common measure changed')
    if Path(cache_path).exists() or not 0 < max_seconds <= 300:
        raise ValueError('A new cache and bounded positive budget required')
    proposal = json.loads((DIRECTORY/PROPOSAL_NAME).read_text(encoding='utf-8'))
    if (proposal['groups'] != GROUPS or proposal['mark'] != str(MARK)
            or proposal['input_marker'] != 1 or not proposal['common_measure_for_both_arities']):
        raise ValueError('Wrong frozen mixed-class scope')
    pins = dict(proposal['reviewed_sources'], **{
        name+'.py': engine_sha, 'hd_bit_syndrome_counted_general.py': COUNTED_SHA})
    for filename, expected in pins.items():
        if sha256((DIRECTORY/filename).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+filename)
    engine = importlib.import_module(name)
    started, stages = perf_counter(), []

    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Extra-arity control exhausted its explicit budget')

    def stage(label):
        row = dict(stage=label, elapsed_seconds=perf_counter()-started)
        stages.append(row)
        if progress:
            progress(row)

    measure = proposal['state_measure']
    if arity == 7:
        prepared = engine.prepare_geometry(check_time=check_time)
        stage('all11440_actual_arity7_fibers')
        answer = engine.full_row_bounds(measure, GROUPS, prepared=prepared, check_time=check_time)
    else:
        answer = engine.full_row_bounds(arity, measure, GROUPS, check_time=check_time)
    stage('extra_arity_ordered_word_contraction')
    old_path = Path(proposal['cache'])
    if not old_path.is_absolute():
        old_path = DIRECTORY/old_path
    if sha256(old_path.read_bytes()).hexdigest() != CACHE_SHA:
        raise ValueError('Frozen mixed8/9 cache changed')
    with np.load(old_path, allow_pickle=False) as old:
        column8, _ = combine_columns(old['quiet8'], old['nonquiet8'], GROUPS, 8, 0.,
                                    lambdas=(F(1, 8), F(7, 8)))
        word8 = log_matrix_power(column8, 16)
        word9 = log_matrix_power(old['column9'], 16)
        norms = {8: combine_row_bounds(np.logaddexp2.reduce(word8[:, 0, :], axis=1), GROUPS)-64*GROUPS*log2(255),
                 9: combine_row_bounds(np.logaddexp2.reduce(word9[:, 0, :], axis=1), GROUPS)-72*GROUPS*log2(255)}
    norms[arity] = answer['combined_word_log2']-8*arity*GROUPS*log2(255)
    count = independent_source_counts()['F_le1']
    prefix = log2(count)-LIMIT*log2(float(MARK))-log2(255)
    mixed = GROUPS*np.logaddexp2.reduce([value/GROUPS for value in norms.values()])
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, column=answer['column_log2_matrices'],
                            rows=np.asarray(answer['rowsets']), words=answer['per_rowset_word_log2'])
    return dict(construction=answer['metadata']['outer_construction'], arity=arity,
        groups=GROUPS, source_class='Fmax<=1', source_count=str(count), mark=str(MARK), input_marker=1,
        per_arity_log2={str(a): prefix+value for a, value in norms.items()},
        normalized_word_log2={str(a): value for a, value in norms.items()},
        common_prefix_log2=prefix, mixed_arities=sorted(norms), mixed_log2=prefix+mixed,
        common_measure_for_all_arities=True, source_count_once=True, projective_normalizations=1,
        metadata=answer['metadata'], costs=answer['costs'], stages=stages,
        elapsed_seconds=perf_counter()-started, reviewed_sources=pins,
        proposal_sha256=PROPOSAL_SHA, cache=str(cache_path),
        cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity', type=int, required=True)
    parser.add_argument('--engine-sha256', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.cache.exists() or args.output.exists():
        parser.error('Do not overwrite a saved control')
    result = run(args.arity, args.engine_sha256, args.cache,
                 progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: result[key] for key in ('per_arity_log2', 'mixed_log2', 'elapsed_seconds')}), flush=True)
