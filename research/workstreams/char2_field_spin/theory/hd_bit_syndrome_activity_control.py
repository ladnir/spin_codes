"""Cached FLOAT test: enforce nine nonzero RS columns before Holder.

For a fixed active row set A, the true column kernel K_A contains the
zero-input term one. Marking every nonzero column by Y>=1 gives
K_Y=1+Y(K_A-1). Every nonzero RS[16,8] row occupies at least nine columns.
Thus each group's word bound acquires Y^-9, outside its Holder root.

If Q=L[B_A^r], B_A>=K_A>=1, then
L[K_Y^r] <= S+Y^r(Q-S), where S=L[1]. The same positive chronological
functional L is used for all row sets and arities. The weaker control
subtracts just the common all-zero-atom component of S. No new regional
provider or numerical fitting runs here. These FLOAT controls are not
directed certificates.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_gate import combine_columns
from hd_post_refresh_quiet_rows import combine_row_bounds, log_matrix_power, moment_atoms

DIRECTORY = Path(__file__).resolve().parent
PROPOSAL = 'hd_bit_syndrome_mixed89_9927_x1_float.json'
PROPOSAL_SHA = 'd94ccf2b9a9bc073b67c90256608fa67ef7be02f44efc47eede301f627c4e2c3'
ARITY7 = 'hd_bit_syndrome_arity7_9927_x1_float.json'
GROUPS = 18


def log_difference(upper, lower):
    """Stable FLOAT log2 difference; reject a negative residual, never clip it."""
    upper, lower = np.broadcast_arrays(np.asarray(upper, dtype=float),
                                       np.asarray(lower, dtype=float))
    if (np.isnan(upper).any() or np.isnan(lower).any()
            or np.isposinf(upper).any() or np.isposinf(lower).any()
            or np.any(upper < lower)):
        raise ArithmeticError('Nonnegative finite-or-zero residual required')
    result = upper.copy()
    finite = np.isfinite(lower)
    with np.errstate(divide='ignore', invalid='raise'):
        result[finite] += np.log2(-np.expm1((lower[finite]-upper[finite])*log(2)))
    return result


def marked_columns(column, baseline, marker, groups=GROUPS):
    marker = F(marker)
    if marker < 1 or type(groups) is not int or groups < 1:
        raise ValueError('Y>=1 and a positive integer Holder exponent required')
    column, baseline = np.asarray(column), np.asarray(baseline)
    if column.shape[-2:] != (2, 2) or baseline.shape != (2, 2):
        raise ValueError('Column batch and common2x2 baseline required')
    if marker == 1:
        return column.copy()
    return np.logaddexp2(baseline, groups*log2(marker)+log_difference(column, baseline))


def constant_columns(measure):
    """The full constant path functional and its common zero-atom subevent."""
    _, empty = moment_atoms(measure)
    full = log_matrix_power(empty, 16)
    zero = np.full((2, 2), -np.inf)
    p = 2.**measure['log2_probability']
    zero[0, 0] = 16*np.log1p(-p)/log(2)
    return full, zero


def norm(column, arity, marker=F(1), groups=GROUPS):
    word = log_matrix_power(column, 16)
    rows = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    return (combine_row_bounds(rows, groups)-8*arity*groups*log2(255)
            -9*groups*log2(marker))


def load_pinned():
    raw = (DIRECTORY/PROPOSAL).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA:
        raise ValueError('Common proposal changed')
    proposal = json.loads(raw)
    a7 = json.loads((DIRECTORY/ARITY7).read_text(encoding='utf-8'))
    if (a7['proposal_sha256'] != PROPOSAL_SHA or a7['arity'] != 7
            or a7['groups'] != GROUPS or a7['mark'] != '9927/10000'
            or a7['source_class'] != 'Fmax<=1' or a7['input_marker'] != 1):
        raise ValueError('Wrong arity7 scope')
    columns, pins = {}, {}
    for receipt in (proposal, a7):
        for name, expected in receipt['reviewed_sources'].items():
            if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
                raise ValueError('Reviewed source changed: '+name)
            pins[name] = expected
        cache = Path(receipt['cache'])
        if not cache.is_absolute():
            cache = DIRECTORY/cache
        if sha256(cache.read_bytes()).hexdigest() != receipt['cache_sha256']:
            raise ValueError('Frozen column cache changed')
        with np.load(cache, allow_pickle=False) as saved:
            if receipt is a7:
                columns[7] = saved['column'].copy()
            else:
                columns[8], _ = combine_columns(saved['quiet8'], saved['nonquiet8'],
                    GROUPS, 8, 0., lambdas=(F(1, 8), F(7, 8)))
                columns[9] = saved['column9'].copy()
    for a, column in columns.items():
        if abs(norm(column, a)-a7['normalized_word_log2'][str(a)]) > 1e-8:
            raise ArithmeticError('Cached baseline reproduction failed')
    return proposal, a7, columns, pins


def run():
    start = perf_counter()
    proposal, a7, columns, pins = load_pinned()
    prefix = a7['common_prefix_log2']
    full, zero = constant_columns(proposal['state_measure'])
    controls = []
    # Mechanism tests, not a free parameter search. a8/a9 remain fixed while
    # the faileda7 class receives the newly enforced nine-column constraint.
    for method, base in (('common_zero_event', zero), ('whole_constant', full)):
        for marker in (F(1), F(2), F(12, 5), F(5, 2)):
            values = {a: norm(marked_columns(column, base, marker), a, marker)
                      if a == 7 else norm(column, a) for a, column in columns.items()}
            mixed = GROUPS*np.logaddexp2.reduce([v/GROUPS for v in values.values()])
            controls.append(dict(method=method, arity7_marker=str(marker),
                arity8_marker='1', arity9_marker='1',
                per_arity_log2={str(a): prefix+v for a, v in values.items()},
                mixed_log2=prefix+mixed))
    return dict(groups=GROUPS, source_class='Fmax<=1', arities=[7, 8, 9],
        mark=proposal['mark'], input_marker=1, active_column_floor_per_group=9,
        common_measure=True, controls=controls, reviewed_sources=pins,
        input_receipts={PROPOSAL: PROPOSAL_SHA,
            ARITY7: sha256((DIRECTORY/ARITY7).read_bytes()).hexdigest()},
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-start, regional_providers=0,
        native_encoder_changed=False, directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a saved control')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('controls', 'elapsed_seconds')}), flush=True)
