"""One full fixed-row-set evaluation of the saved 99/100 output mark.

The prescribed Young weights are unchanged. This is still only the complete
r18/arity8/Fmax1 class, in FLOAT. No new regional provider runs here.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_post_refresh_full_rows import full_row_bounds, LAMBDAS
from hd_post_refresh_moment_gate import check_source_pins
from hd_two_sum_closed_kernel import branch_counts


DIRECTORY = Path(__file__).resolve().parent
PROBE = 'hd_post_refresh_mark_probe_float.json'
PROBE_SHA = '7495e07911e61ae54321d743d0dbcea87dfb67c903bb959540b07b5b506568c3'
SOURCE_PINS = {
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
    'hd_post_refresh_full_rows.py': '982c9a1abd836939d6e822d03d969d29961de9df8e81c86a15c23f23834445dc',
}


def run(log2_x, cache_path):
    if not np.isfinite(log2_x) or log2_x < 0:
        raise ValueError('Finite scalar X>=1 required')
    if Path(cache_path).exists():
        raise ValueError('Do not overwrite a saved branch bank')
    check_source_pins()
    for name, pin in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != pin:
            raise ValueError('Reviewed source changed: '+name)
    raw = (DIRECTORY/PROBE).read_bytes()
    if sha256(raw).hexdigest() != PROBE_SHA:
        raise ValueError('Output-mark probe changed')
    probe = json.loads(raw)
    data, = [row for row in probe['controls'] if row['mark'] == '99/100']
    if not data['complete']:
        raise ValueError('Complete output-mark control required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started > 120:
            raise TimeoutError('Full row-set gate exceeded two minutes')
    def progress(stage):
        print(json.dumps(dict(stage=stage, elapsed_seconds=perf_counter()-started)), flush=True)
    counts = [branch_counts(8, b, 16) for b in range(9)]
    logs = np.array([[log2(row[i]) if row[i] else -np.inf for row in counts] for i in (0, 2, 1)])
    answer = full_row_bounds(logs, 18, data['state_measure'], log2_x,
        batch_size=128, check_time=check_time, progress=progress)
    # These two branch matrices are independent of the activity marker X.
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, full=answer['branch_column_log2_matrices']['full'],
            holes=answer['branch_column_log2_matrices']['holes'],
            rowsets=np.asarray(answer['rowsets']), output_mark='99/100',
            probe_sha256=PROBE_SHA, source_sha256=SOURCE_PINS['hd_post_refresh_full_rows.py'])
    prefix = data['prefix_log2']-1296*log2_x
    quiet_word = __import__('hd_post_refresh_quiet_rows').log_matrix_power(
        answer['branch_column_log2_matrices']['quiet'], 16)
    quiet_values = np.logaddexp2.reduce(quiet_word[:, 0, :], axis=1)
    quiet_log = float(18*np.logaddexp2.reduce(quiet_values/18))
    return dict(construction=probe['construction'], scope=probe['scope'],
        mark='99/100', activity_marker_log2=log2_x,
        young_weights=list(map(str, LAMBDAS)), probe_artifact=PROBE, probe_sha256=PROBE_SHA,
        reviewed_source_sha256=SOURCE_PINS, branch_cache=str(cache_path),
        branch_cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        prefix_log2=prefix, quiet_component_log2=prefix+quiet_log,
        combined_word_log2=answer['combined_word_log2'],
        full_class_log2=prefix+answer['combined_word_log2'],
        fixed_rowsets=True, rowset_count=len(answer['rowsets']),
        raw_regional_banks_built=0, directed=False, whole_code_certificate=False,
        complete=True, elapsed_seconds=perf_counter()-started,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log2-x', required=True, type=float)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Do not overwrite a frozen result or branch cache')
    result = run(args.log2_x, args.cache)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, allow_nan=False))
