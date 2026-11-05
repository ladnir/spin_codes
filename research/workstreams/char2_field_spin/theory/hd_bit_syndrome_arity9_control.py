"""One bounded FLOAT a9/a8+a9 coverage gate using the frozen common measure.

This control changes no encoder and computes no fresh regional provider.
The saved FLOAT proposal supplies the same common measure for both arities.
The a9 engine includes the actual zero-outside fiber; a8 keeps its reviewed
two-branch Young bound. A pass is numerical evidence for this mixed class,
not a directed replay or a whole-code certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_arity9_gate import full_row_bounds, prepare_geometry
from hd_bit_syndrome_directed_recheck import independent_source_counts
from hd_bit_syndrome_gate import combine_columns
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power

DIRECTORY = Path(__file__).resolve().parent
PROPOSAL = DIRECTORY/'hd_bit_syndrome_orbit_993_x0_float.json'
PROPOSAL_SHA = '4d5d2374edca033663a18de558b1eef56751e9e4fc6b8081df5ec13c0f89e369'
BRANCH_SHA = '5d68a6899c5ba9ae88836a5347dc78f5a2ab10914286a6337ceda933c745531c'
MARK, GROUPS, LIMIT = F(993, 1000), 18, 734003


def run(engine_sha256, *, output_cache, max_seconds=180, progress=None):
    engine = DIRECTORY/'hd_bit_syndrome_arity9_gate.py'
    if len(engine_sha256) != 64 or sha256(engine.read_bytes()).hexdigest() != engine_sha256:
        raise ValueError('Reviewed arity9 engine changed')
    if sha256(PROPOSAL.read_bytes()).hexdigest() != PROPOSAL_SHA:
        raise ValueError('Frozen common-measure proposal changed')
    if not 0 < max_seconds <= 300 or Path(output_cache).exists():
        raise ValueError('Bounded budget and a new output cache required')
    proposal = json.loads(PROPOSAL.read_text(encoding='utf-8'))
    if (proposal['mark'] != str(MARK) or proposal['active_groups'] != GROUPS
            or proposal['log2_activity_marker'] != 0 or proposal['whole_code_certificate']):
        raise ValueError('Wrong common-measure proposal scope')
    measure = proposal['state_measure']
    cache_path = Path(proposal['branch_cache'])
    if not cache_path.is_absolute():
        cache_path = DIRECTORY/cache_path
    if sha256(cache_path.read_bytes()).hexdigest() != BRANCH_SHA:
        raise ValueError('Frozen arity8 branches changed')
    started, stages = perf_counter(), []

    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Arity9 bounded gate exhausted its budget')

    def stage(name):
        record = dict(stage=name, elapsed_seconds=perf_counter()-started)
        stages.append(record)
        if progress:
            progress(record)

    prepared = prepare_geometry(check_time=check_time)
    stage('all11440_actual_fibers_prepared')
    answer = full_row_bounds(measure, GROUPS, prepared=prepared, check_time=check_time)
    stage('ordered_arity9_word_contraction')
    source_counts = independent_source_counts()
    common = log2(source_counts['F_le1']) - LIMIT*log2(float(MARK)) - log2(255)
    normalized9 = answer['combined_word_log2'] - 72*GROUPS*log2(255)
    with np.load(cache_path, allow_pickle=False) as old:
        combined8, _ = combine_columns(old['quiet'], old['nonquiet'], GROUPS, 8, 0.,
                                       lambdas=(F(1, 16), F(15, 16)))
    words8 = log_matrix_power(combined8, 16)
    normalized8 = combine_row_bounds(np.logaddexp2.reduce(words8[:, 0, :], axis=1), GROUPS)-64*GROUPS*log2(255)
    mixed = GROUPS*np.logaddexp2(normalized8/GROUPS, normalized9/GROUPS)
    with Path(output_cache).open('xb') as handle:
        np.savez_compressed(handle, column=answer['column_log2_matrices'],
            rowsets=np.asarray(answer['rowsets']), words=answer['per_rowset_word_log2'])
    return dict(construction=answer['metadata']['outer_construction'], groups=GROUPS,
        source_class='Fmax<=1', arity9_log2=common+normalized9, arity8_log2=common+normalized8,
        mixed_arity8_or9_log2=common+mixed, mark=str(MARK), input_marker=1,
        all_arity9_zero_outside_fibers=True, common_measure_for_both_arities=True,
        mixed_method='18th power of the sum of normalized per-arity18th roots',
        source_count_once=True, projective_normalizations=1, metadata=answer['metadata'],
        costs=answer['costs'], stages=stages, elapsed_seconds=perf_counter()-started,
        a9_engine_sha256=engine_sha256, proposal_sha256=PROPOSAL_SHA,
        a8_branch_cache_sha256=BRANCH_SHA, output_cache=str(output_cache),
        output_cache_sha256=sha256(Path(output_cache).read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine-sha256', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.cache.exists() or args.output.exists():
        parser.error('Do not overwrite a saved control')
    answer = run(args.engine_sha256, output_cache=args.cache,
                 progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: answer[k] for k in ('arity9_log2', 'arity8_log2',
        'mixed_arity8_or9_log2', 'elapsed_seconds')}), flush=True)
