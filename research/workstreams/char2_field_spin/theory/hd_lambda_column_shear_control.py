"""Selected r18/Fmax<=1 FLOAT gate for the deterministic column shear.

Fix one original row set A of size a=1 or2 throughout all16 columns.
In column j, use G_j=Q_j H_lambda P_j[:,A], with lambda=0xbc, inverse
input order, and output xor-j order. For each even j, row parity selects
which endpoint receives the other: (u,v)->(u+v,v) or (u,u+v).
The final independent nonzero diagonal is AFTER this shear.

The32 output forms of a pair act on2a original input variables. Their
GF4 relation kernels give exact raw counts n_k, original input weights
w_k, and BOTH output masks jointly. We do not multiply independent
column support laws or use the earlier fresh-ratio pair law.

Under the fixed positive regional matrix measure, J_k is the r-th
moment of the resulting32-coordinate support monomial. Entrywise,
P(X)=sum_k n_k J_k**(1/r) X**w_k; coefficientwise Holder bounds the
pair's r-group moment by P(X)**r. Eight ordered matrix-polynomial
products retain degrees >=9*a*r. Finally sum fixed-A roots, divide
each group by255**(8*a), and apply the common source/Chernoff/projective
prefix ONCE. This is a positive FLOAT calculation, not a certificate.

The RS comparison assumes ideal independent D1 coefficients before the
whole pair map. One common pair scale commutes through the shear and
can be absorbed into its final D2. Conditional on that scale, the key
change is a product-uniform bijection. Two separate column gauges do
NOT commute through this shear. Compared with the old normalized
implementation, the mate's pivot must be restored: eight additional
D1 products per superblock. This file performs no native operations.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_lambda_column_shear_ranks import paired_rows, SCHEDULE
from hd_quaternary_outer_relations import category_data, column_masks
from hd_pair_identity_kernel import mask_log_matrices, branch_log_matrices
from hd_pair_scheduled_control import load_measure
from hd_pair_state_tail import word_tail

DIRECTORY = Path(__file__).resolve().parent
INPUT = 'hd_mds_exact_empty_mixed789_9927_pinned_float.json'
INPUT_SHA = 'c89ef4973de29007fe293817d2e3b9f582246a45d19ca19c7fcce2554e7d33b8'
GROUPS = 18
PINS = {
    'hd_lambda_column_shear_ranks.py': '4067013201cadef412ae647cc81e394e90b49d9403dfdd4713a739a191a96fa8',
    'hd_quaternary_outer_relations.py': 'd8be4e7a1ab77b98c76695f13c1bc9a9cf867df2aae6408ac3e08c6fc12346c0',
    'hd_pair_identity_kernel.py': '93f3a9053a72bed39db6b00a32d2f2b7b060a35ce19ce814dd856d096a796961',
    'hd_pair_state_tail.py': '2af63d37d228a76a26853d7550815de5f6c0bc5af05908a5b009e14aa9175858',
    'hd_pair_scheduled_control.py': 'cd0ee6b03e794785dee26502940b3f68395cd9a1fd6f27dc6f3470a490730dd4',
    'hd_pair_scalar_control.py': '088f33f37dc70b9da7d4649d7f56e8657fd84d4945c909bf49235cc20d5c3f61',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
    'hd_inverse_column_control.py': '40668f4bf673ef0859be429a82119769df169408d74aaaa52dfbb68787b9716b',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
    'hd_pair_state_kernel.py': 'b28142e19e51a906ae742128325d1e7557dfb425624c6da45778a2a75753e195',
}


def dependencies():
    """Pin the cached bank, its original providers, and reused code."""
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Frozen common-measure receipt changed')
    saved = json.loads(raw)
    pins = dict(saved['reviewed_sources'])
    for name, digest in PINS.items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting frozen provenance: '+name)
        pins[name] = digest
    for name, digest in pins.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    pins[INPUT] = INPUT_SHA
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return pins


def _arity(arity):
    if type(arity) is not int or arity not in (1, 2):
        raise ValueError('This gate supports exactly arities1 and2')


def pair_data(source, column, *, q=256):
    """Exact joint categories; bit0..15 is first, bit16..31 is mate."""
    source = tuple(source)
    _arity(len(source))
    if type(column) is not int or column not in range(0, 16, 2):
        raise ValueError('An even physical column0..14 is required')
    left, right = paired_rows(source, column, 1, schedule=SCHEDULE)
    categories = category_data(2*len(source), q=q)
    masks = column_masks(left+right, 2*len(source))
    counts = categories['category_counts']
    weights = categories['category_input_weights']
    if sum(counts) != q**(2*len(source)):
        raise ArithmeticError('Joint raw input mass changed')
    for weight, mask, count in zip(weights, masks, counts):
        if count and (weight == 0) != (mask == 0):
            raise ArithmeticError('The shear must preserve the unique zero input')
    return dict(source=source, column=column, mate=column+1,
        first=tuple(mask & 65535 for mask in masks),
        second=tuple(mask >> 16 for mask in masks),
        input_weights=weights, raw_counts=counts,
        input_masks=categories['category_input_masks'])


def branch_data(arity, sources=None, *, check_time=None):
    """Arrays(8,number_of_A,categories); source labels stay fixed.

    Supplying a subset is useful for bounded batches. Only run() decides
    whether all16 or120 labelled row sets have actually completed.
    """
    _arity(arity)
    sources = tuple(combinations(range(16), arity)) if sources is None else tuple(map(tuple, sources))
    if (not sources or len(set(sources)) != len(sources)
            or any(len(source) != arity or tuple(sorted(source)) != source
                   or len(set(source)) != arity
                   or any(type(i) is not int or not 0 <= i < 16 for i in source)
                   for source in sources)):
        raise ValueError('Distinct increasing original row sets of the requested arity required')
    categories = category_data(2*arity)
    counts, weights = categories['category_counts'], categories['category_input_weights']
    branch_count = len(counts)
    first = np.empty((8, len(sources), branch_count), dtype=np.uint16)
    second = np.empty_like(first)
    for block in range(8):
        if check_time is not None:
            check_time()
        for index, source in enumerate(sources):
            data = pair_data(source, 2*block)
            if data['raw_counts'] != counts or data['input_weights'] != weights:
                raise ArithmeticError('A schedule changed the joint category partition')
            first[block, index] = data['first']
            second[block, index] = data['second']
    logs = np.array([log2(n) if n else -np.inf for n in counts])
    return dict(sources=sources, branch_count=branch_count, first=first, second=second,
        raw_counts=counts, log_counts=np.broadcast_to(logs, (len(sources), branch_count)),
        input_weights=np.broadcast_to(np.array(weights, dtype=np.int8), (len(sources), branch_count)))


def normalized_roots(tails, arity):
    _arity(arity)
    values = np.asarray(tails, dtype=float)
    if values.ndim != 1 or np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('One finite or negative-infinite log tail per labelled A required')
    return values/GROUPS-8*arity*log2(255)


def root_sum_bound(roots, prefix):
    roots = np.asarray(roots, dtype=float)
    if (roots.ndim != 1 or not len(roots) or np.any(np.isnan(roots))
            or np.any(np.isposinf(roots)) or not isfinite(prefix)):
        raise ValueError('Nonempty valid root array and finite common prefix required')
    return prefix+GROUPS*float(np.logaddexp2.reduce(roots))


def _json_log(value):
    """JSON null denotes exact zero, whose log is negative infinity."""
    return None if np.isneginf(value) else float(value)


def run(*, arities=(1, 2), batch_size=8, max_seconds=120, progress=None):
    """Root-owned bounded selected gate; no provider or optimizer calls."""
    arities = tuple(arities)
    if (not arities or len(set(arities)) != len(arities)
            or any(type(a) is not int or a not in (1, 2) for a in arities)
            or type(batch_size) is not int or not 1 <= batch_size <= 15
            or not isfinite(max_seconds) or not 0 < max_seconds <= 600):
        raise ValueError('Distinct arities1/2, batch1..15, and time budget(0,600] required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('Column-shear selected gate exceeded its budget')
    pins = dependencies()
    saved, measure, prefix = load_measure()
    if (saved['groups'] != GROUPS or saved['source_class'] != 'Fmax<=1'
            or saved['mark'] != '9927/10000'):
        raise ValueError('Wrong regional source scope')
    rows, controls, all_roots = [], [], []
    complete, timeout = True, None
    try:
        check_time()
        table = mask_log_matrices(measure, GROUPS)
        for arity in arities:
            sources = tuple(combinations(range(16), arity))
            roots = []
            for begin in range(0, len(sources), batch_size):
                check_time()
                selected = sources[begin:begin+batch_size]
                data = branch_data(arity, selected, check_time=check_time)
                blocks = []
                for block in range(8):
                    check_time()
                    matrices = branch_log_matrices(data['first'][block].reshape(-1),
                        data['second'][block].reshape(-1), table)
                    blocks.append(matrices.reshape(len(selected), data['branch_count'], 2, 2))
                tails = word_tail(blocks, data['log_counts'], data['input_weights'],
                                  GROUPS, 9*arity*GROUPS, check_time=check_time)
                values = normalized_roots(tails, arity)
                roots.extend(values.tolist())
                rows.extend(dict(arity=arity, source=list(source), tail_log2=_json_log(tail),
                    normalized_group_root_log2=_json_log(root))
                    for source, tail, root in zip(selected, tails, values))
                if progress:
                    progress(dict(stage='completed_rowsets', arity=arity, completed=len(roots),
                        requested=len(sources), elapsed_seconds=perf_counter()-start))
            bound = root_sum_bound(roots, prefix)
            controls.append(dict(arity=arity, labelled_rowsets=len(sources),
                expected_rowsets=comb(16, arity), complete=True, floor=9*arity*GROUPS,
                maximum_degree=16*arity*GROUPS, branch_count=data['branch_count'],
                exact_raw_mass_per_pair=str(256**(2*arity)), bound_log2=_json_log(bound),
                largest_fixed_A=list(sources[int(np.argmax(roots))])))
            all_roots.extend(roots)
            if progress:
                progress(dict(stage='completed_arity', **controls[-1]))
    except TimeoutError as error:
        complete, timeout = False, str(error)
    if dependencies() != pins:
        raise ValueError('Frozen sources or bank changed during this calculation')
    mixed = root_sum_bound(all_roots, prefix) if complete else None
    return dict(construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        arithmetic='positive-FLOAT-log2-not-directed', complete=complete, timeout=timeout,
        requested_arities=list(arities), completed_arities=[row['arity'] for row in controls],
        all_labelled_arities=[row['arity'] for row in controls], groups=GROUPS,
        source_class=saved['source_class'], source_count=saved['source_count'], mark=saved['mark'],
        common_prefix_log2=prefix, source_and_projective_prefix_applied_once=True,
        schedule=SCHEDULE, normal=1, column_pairs=[[j,j+1] for j in range(0,16,2)],
        original_source_A_fixed=True, exact_joint_GF4_categories=True,
        independent_column_support_assumption=False, fresh_ratio_law_used=False,
        original_input_weight_floor=True, separate_column_activity_floor_retained=False,
        full_independent_D1_comparison=True, permitted_gauges_per_column_pair=1,
        restored_D1_products_per_superblock=8, final_D2_after_shear=True,
        common_measure=measure, controls=controls, fixed_A_rows=rows,
        requested_arity_union_log2=_json_log(mixed) if mixed is not None else None,
        partial_rows_are_not_a_complete_class_bound=not complete,
        input_receipt=INPUT, input_sha256=INPUT_SHA, source_pins=pins,
        provider_runs=0, optimizations=0, directed=False, whole_code_certificate=False,
        native_implementation_changed=False, elapsed_seconds=perf_counter()-start)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arities', nargs='+', type=int, choices=(1,2), default=(1,2))
    parser.add_argument('--batch-size', type=int, default=8)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(arities=args.arities, batch_size=args.batch_size, max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=result['complete'], controls=result['controls'],
        requested_arity_union_log2=result['requested_arity_union_log2'],
        elapsed_seconds=result['elapsed_seconds']), allow_nan=False), flush=True)
