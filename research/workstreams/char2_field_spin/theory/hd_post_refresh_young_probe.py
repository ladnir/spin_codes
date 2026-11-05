"""Saved-mark full-row comparison and four cheap Young-weight witnesses.

One shared set of branch matrices is built. Changing Young weights does
not rebuild any regional bank or coefficient dynamic program.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_post_refresh_full_mark_gate import DIRECTORY, PROBE, PROBE_SHA, SOURCE_PINS
from hd_post_refresh_full_rows import full_row_bounds, combine_columns
from hd_post_refresh_quiet_rows import log_matrix_power, combine_row_bounds
from hd_post_refresh_moment_gate import check_source_pins
from hd_two_sum_closed_kernel import branch_counts


YOUNG_WEIGHTS = ((F(1, 3),)*3, (F(1, 5), F(3, 5), F(1, 5)),
                 (F(1, 10), F(4, 5), F(1, 10)), (F(2, 5), F(1, 5), F(2, 5)))


def run(mark, log2_x, cache_path):
    check_source_pins()
    if not np.isfinite(log2_x) or log2_x < 0 or Path(cache_path).exists():
        raise ValueError('Finite X>=1 and a new cache path required')
    for name, pin in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != pin:
            raise ValueError('Reviewed source changed: '+name)
    raw = (DIRECTORY/PROBE).read_bytes()
    if sha256(raw).hexdigest() != PROBE_SHA:
        raise ValueError('Pinned mark probe changed')
    probe = json.loads(raw)
    data, = [row for row in probe['controls'] if row['mark'] == mark and row['complete']]
    started = perf_counter()
    def check_time():
        if perf_counter()-started > 120:
            raise TimeoutError('Full row-set witness probe exceeded two minutes')
    def progress(stage):
        print(json.dumps(dict(stage=stage, elapsed_seconds=perf_counter()-started)), flush=True)
    counts = [branch_counts(8, b, 16) for b in range(9)]
    logs = np.array([[log2(row[i]) if row[i] else -np.inf for row in counts] for i in (0, 2, 1)])
    answer = full_row_bounds(logs, 18, data['state_measure'], log2_x,
        batch_size=128, check_time=check_time, progress=progress)
    branches = answer['branch_column_log2_matrices']
    prefix = data['prefix_log2']-1296*log2_x
    def evaluate(matrix):
        word = log_matrix_power(matrix, 16)
        return prefix+combine_row_bounds(np.logaddexp2.reduce(word[:, 0, :], axis=1), 18)
    results = []
    for weights in YOUNG_WEIGHTS:
        check_time()
        matrix, _ = combine_columns(branches['quiet'], branches['full'], branches['holes'],
            18, 8, 8, log2_x, lambdas=weights)
        row = dict(young_weights=list(map(str, weights)), log2_screen=evaluate(matrix))
        results.append(row)
        print(json.dumps(row), flush=True)
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, **branches, rowsets=np.asarray(answer['rowsets']),
            output_mark=mark, activity_marker_log2=log2_x, probe_sha256=PROBE_SHA)
    return dict(construction=probe['construction'], scope=probe['scope'],
        output_mark=mark, activity_marker_log2=log2_x, prefix_log2=prefix,
        quiet_component_log2=evaluate(branches['quiet']), controls=results,
        raw_regional_banks_built=0, branch_banks_built=1, young_witness_count=len(results),
        directed=False, whole_code_certificate=False, complete=True,
        probe_sha256=PROBE_SHA, reviewed_source_sha256=SOURCE_PINS,
        cache_artifact=str(cache_path), cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mark', required=True)
    parser.add_argument('--log2-x', required=True, type=float)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Do not overwrite a frozen result or cache')
    result = run(args.mark, args.log2_x, args.cache)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, allow_nan=False))
