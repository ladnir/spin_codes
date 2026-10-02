"""Fresh sparse-occupancy search with separately authenticated improved counts.

This is a research receipt, not a complete certificate. Each successful row
has an exhaustive support-box cover and fresh outward verification. Missing
or failed occupancies remain explicitly outside the claim.
"""
import argparse
from fractions import Fraction as Q
import os
from pathlib import Path

from flint import arb, ctx
import sparse_hill_q1 as q1


def validate(occupancies, distance, tilts, precision, target_bits, max_splits, output):
    if (not occupancies or len(set(occupancies)) != len(occupancies)
            or any(type(q) is not int or not 2 <= q <= 32 for q in occupancies)
            or not 0 < Q(distance) < Q(1, 2)
            or not tilts or len(set(map(Q, tilts))) != len(tilts) or any(Q(t) <= 0 for t in tilts)
            or type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 1
            or type(max_splits) is not int or max_splits < 0
            or output.exists()):
        raise ValueError('distinct q=2..32, positive tilts, valid limits, and fresh output required')


def run(occupancies, distance, tilts, precision, target_bits, max_splits, output,
        incidence=False):
    validate(occupancies, distance, tilts, precision, target_bits, max_splits, output)
    sparse = q1._sparse_interface()
    if incidence:
        from outer_hill_incidence import authenticated_bch_cdf
        counts, premises = authenticated_bch_cdf()
    else:
        counts, premises = q1.authenticated_counts(sparse, True)
    counts, floor = sparse.checked_counts(counts)
    threshold = (Q(distance)*(1 << 21)).__floor__()
    record = dict(schema='packed-gl32-sparse-hill-1', K=1 << 20, N=1 << 21,
        group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
        mixing='independent uniform GL32 per group/block; no additional GF16 stage',
        routing='independent shared column shuffle per group and independent regional shuffles',
        inner=dict(t=128, s=19, updates=2), distance=str(Q(distance)), threshold=threshold,
        precision=precision, target_bits=target_bits, tilts=list(tilts), max_splits=max_splits,
        joint_return_through=2, lazy_density_through=4,
        count_sha256=sparse.count_hash(counts), count_premises=premises, support_min=floor,
        requested=list(occupancies), results=[], process_id=os.getpid(),
        count_refinement='H5 exact + incidence' if incidence else 'H5 exact',
        complete_requested=False, complete_sparse_prefix=False,
        note='Only successful listed occupancies; q1 and dense occupancies are outside this receipt.')
    sparse.save(output, record)
    args = sparse.old_args.build_args(max(occupancies), tilts, precision, max_splits, target_bits, None, 2)
    args.exact_feedback = True
    args.analytic_gradient = True
    args.joint_return_through = 2
    args.lazy_density_through = 4
    operators = sparse.occupancy_birth_classes.build_operators(args)
    terminal = sparse.np.ones(next(iter(operators.values()))[0][0].nrows())
    from support_cover import cover
    total = arb(0)
    for q in occupancies:
        args.groups = q
        ctx.prec = precision
        upper, details = cover(args, operators, counts, terminal, cutoff=threshold, support_min=floor)
        passed = upper is not None and upper < arb(2)**(-target_bits)
        record['results'].append(dict(occupancy=q, passed=passed,
            upper=sparse.endpoint(upper) if upper is not None else None, details=details))
        if passed:
            total = sparse.up(total+upper)
            print('SPARSE HILL q', q, 'margin', -upper.log()/arb(2).log(), flush=True)
        else:
            print('SPARSE HILL q', q, 'not closed within split budget', flush=True)
        record['successful_aggregate_upper'] = sparse.endpoint(total) if total > 0 else None
        record['complete_requested'] = (len(record['results']) == len(occupancies)
            and all(row['passed'] for row in record['results']))
        sparse.save(output, record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies', type=int, nargs='+', default=list(range(2, 33)))
    parser.add_argument('--distance', default='.10')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=46)
    parser.add_argument('--max-splits', type=int, default=64)
    parser.add_argument('--tilts', nargs='+', default=['.00032', '.00064', '.001', '.0016', '.0032',
        '.0064', '.008', '.012', '.016', '.024', '.032'])
    parser.add_argument('--incidence', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print('SPARSE HILL PID', os.getpid(), flush=True)
    run(args.occupancies, args.distance, args.tilts, args.precision,
        args.target_bits, args.max_splits, args.output, args.incidence)


if __name__ == '__main__':
    main()
