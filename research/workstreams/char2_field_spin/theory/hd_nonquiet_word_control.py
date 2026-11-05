"""Cached FLOAT control of a positive first-nonquiet word decomposition.

Outside a separately bounded bad-outer event, each source group must have
a nonquiet column. Expand at its first such column BEFORE applying the
group Holder exponent. Convexity gives16^(r-1) times an ordered sum of
quiet-prefix, nonquiet-column, full-suffix moment matrices. Never subtract
one upper matrix from another. This control covers only r18/a8/F<=1;
the inverse-scheduled ensemble needs its own global setup-event bound.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path

import numpy as np

from hd_inverse_column_control import INPUT, INPUT_SHA, schedule
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matmul

DIRECTORY = Path(__file__).resolve().parent


def first_nonquiet_words(quiet, nonquiet, full, rows, permutations):
    rows = tuple(tuple(int(i) for i in row) for row in rows)
    index = {row: i for i, row in enumerate(rows)}
    if quiet.shape != nonquiet.shape or quiet.shape != full.shape or quiet.shape != (len(rows), 2, 2):
        raise ValueError('Matched per-row-set2x2 matrices required')
    before = np.full(quiet.shape, -np.inf)
    before[:, 0, 0] = before[:, 1, 1] = 0.
    after = np.full(quiet.shape, -np.inf)
    for p in permutations:
        if sorted(p) != list(range(16)):
            raise ValueError('Sixteen-coordinate permutation required')
        selected = [index[tuple(sorted(p[i] for i in row))] for row in rows]
        after = np.logaddexp2(log_matmul(after, full[selected]),
                             log_matmul(before, nonquiet[selected]))
        before = log_matmul(before, quiet[selected])
    return np.logaddexp2.reduce(after[:, 0, :], axis=1)


def run():
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Frozen input changed')
    saved = json.loads(raw)
    cache = Path(saved['cache'])
    if not cache.is_absolute():
        cache = DIRECTORY/cache
    if sha256(cache.read_bytes()).hexdigest() != saved['cache_sha256']:
        raise ValueError('Column cache changed')
    r, a = saved['groups'], 8
    results = []
    with np.load(cache, allow_pickle=False) as data:
        for name, maps in (('static', (tuple(range(16)),)*16), ('inverse', schedule())):
            word = first_nonquiet_words(data['quiet8'], data['nonquiet8'],
                                       data['column8'], data['rowsets8'], maps)
            log_bound = (saved['common_prefix_log2'] + combine_row_bounds(word, r)
                         +4*(r-1) -8*a*r*log2(255))
            results.append(dict(schedule=name, restricted_first_moment_log2=log_bound))
    return dict(groups=r, arity=a, source_class=saved['source_class'], mark=saved['mark'],
        nonquiet_columns_per_group=1, decomposition='first nonquiet column before group power',
        convexity_cost_log2=4*(r-1), results=results,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        input_receipt=INPUT, input_sha256=INPUT_SHA,
        cache_sha256=saved['cache_sha256'],
        additional_bad_setup_probability_required=True, bad_setup_charge_included=False,
        static_arity8_bad_setup_bound_available=False, directed=False,
        whole_code_certificate=False, native_changed=False, regional_providers_run=0)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite saved evidence')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)
