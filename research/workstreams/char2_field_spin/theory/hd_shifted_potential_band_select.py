"""Select a complete subinterval of an existing FLOAT screen without recomputation."""
from copy import deepcopy
from math import isfinite
from pathlib import Path

import numpy as np

from hd_shifted_swap_screen import CONSTRUCTION


def select_interval(screen, first, last, *, threshold=-80.):
    """Preserve each recorded weight and exact witness; do not certify saved values."""
    if (screen.get('construction') != CONSTRUCTION or screen.get('directed') is not False
            or screen.get('every_weight_evaluated') is not True or not isfinite(threshold)
            or not isinstance(first, int) or not isinstance(last, int)
            or not screen['requested_weight_range'][0] <= first <= last <= screen['requested_weight_range'][1]):
        raise ValueError('A complete FLOAT screen and a contained integer interval are required')
    rows = [row for row in screen['results'] if first <= row['active_groups'] <= last]
    if [row['active_groups'] for row in rows] != list(range(first, last+1)):
        raise ValueError('Every selected integer weight must appear exactly once in order')
    if any(not isfinite(row['log2_upper_screen']) for row in rows):
        raise ValueError('Finite per-weight FLOAT screen values are required')
    result = deepcopy(screen)
    union = float(np.logaddexp2.reduce([row['log2_upper_screen'] for row in rows]))
    failed = [row['active_groups'] for row in rows if row['log2_upper_screen'] > threshold]
    result.update(requested_weight_range=[first, last], weights_evaluated=len(rows),
                  results=deepcopy(rows), failed_weights=failed, threshold=threshold,
                  all_weights_pass_threshold=not failed, union_log2_upper_screen=union,
                  union_passes_threshold=union <= threshold, directed=False, whole_code_certificate=False,
                  selected_from_weight_range=screen['requested_weight_range'],
                  new_fft_tables_built=0, new_witness_searches=0)
    return result


if __name__ == '__main__':
    import argparse, json
    from hashlib import sha256
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--threshold', type=float, default=-80.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    raw = args.source.read_bytes()
    result = select_interval(json.loads(raw), args.first, args.last, threshold=args.threshold)
    result['source_screen_sha256'] = sha256(raw).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, separators=(',', ':'), allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), requested_weight_range=result['requested_weight_range'],
                          weights_evaluated=result['weights_evaluated'],
                          union_log2_upper_screen=result['union_log2_upper_screen'],
                          union_passes_threshold=result['union_passes_threshold'])), flush=True)
