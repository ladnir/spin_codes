"""Fresh partial distance bounds for independently shuffled row pairs.

Every run regenerates the expected outer CDF and conditional inner bounds.
Only the listed occupancies are covered; no older certificate is imported.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from flint import arb, ctx

import pairwise_support
import sparse_cover
import occupancy_birth_classes
import single_group
import occupancy_rank


def run(occupancies, tilts, precision=256, max_splits=64, target_bits=48,
        threshold=209715, updates=4, refined_counts=False, coupled_counts=False,
        output=None,shell_cdf=False):
    if (not occupancies or len(set(occupancies)) != len(occupancies)
            or any(type(q) is not int or not 1 <= q <= 2048 for q in occupancies)
            or not tilts or any(Q(t) <= 0 for t in tilts)
            or type(precision) is not int or precision < 128
            or type(max_splits) is not int or max_splits < 0
            or type(target_bits) is not int or target_bits < 40
            or type(threshold) is not int or not 0 <= threshold < 1 << 21
            or type(updates) is not int or updates not in (2,3,4)):
        raise ValueError('valid distinct occupancies and proof/search parameters required')
    ctx.prec = precision
    if shell_cdf and not (refined_counts and coupled_counts):
        raise ValueError('shell CDF refinement requires both count-refinement options')
    counts, _ = (pairwise_support.shell_refined_counts() if shell_cdf else
                 pairwise_support.pairwise_counts(refined_counts,coupled_counts))
    args = sparse_cover.build_args(max(occupancies), tilts, precision,
                                   max_splits, target_bits, None, updates)
    args.exact_feedback = True
    args.analytic_gradient = True
    args.proposal_buffer_bits = 2
    args.joint_return_through = 3
    args.lazy_density_through = 6
    operators = occupancy_birth_classes.build_operators(args)
    terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
    results = []
    for q in occupancies:
        print('PAIRWISE GF16 COVER updates', updates, 'occupancy', q,
              'cutoff', threshold, flush=True)
        if q == 1:
            best = [arb(1)]*257
            for tilt in tilts:
                exact = operators[tilt,'1'][0]
                moments = single_group.support_moments(exact[0], exact[1])
                multiplier = (occupancy_rank.aq(Q(tilt))*threshold).exp()
                best = [min(old,occupancy_rank.up(value*multiplier))
                        for old,value in zip(best,moments)]
            upper = occupancy_rank.up(2048*single_group.fold_cdf(counts,best))
        else:
            args.groups = q
            upper = sparse_cover.sparse.cover(args, operators, {'1': counts}, terminal,
                                               cutoff=threshold)
        if upper is not None:
            print('PAIRWISE VERIFIED occupancy', q, 'margin',
                  -upper.log()/arb(2).log(), flush=True)
        results.append(dict(occupancy=q, upper=None if upper is None
                            else [int(v) for v in upper.upper().man_exp()]))
        if output:
            record = dict(schema='pairwise-gf16-sparse-screen-1',
                          ensemble='pairwise4-gf16-r'+str(updates),
                          precision=precision, updates=updates,
                          occupancies=occupancies, threshold=threshold, tilts=tilts,
                          max_splits=max_splits, target_bits=target_bits,
                          refined_counts=refined_counts, coupled_counts=coupled_counts,
                          shell_cdf=shell_cdf,
                          joint_return_through=3, lazy_density_through=6,
                          results=results,
                          note='Expected outer CDF over two independent pair shuffles per group. '
                               'Listed occupancies only; not a full-code certificate.')
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(record, indent=2)+'\n')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies', type=int, nargs='+', default=[1,23,32,48])
    parser.add_argument('--tilts', nargs='+',
                        default=['.00024','.00028','.00032','.0004','.001','.0032','.008','.016','.032','.064','.096'])
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--max-splits', type=int, default=64)
    parser.add_argument('--target-bits', type=int, default=48)
    parser.add_argument('--threshold', type=int, default=209715)
    parser.add_argument('--updates', type=int, choices=(2,3,4), default=4)
    parser.add_argument('--refined-counts', action='store_true')
    parser.add_argument('--coupled-counts', action='store_true')
    parser.add_argument('--shell-cdf',action='store_true',help='Intersect pair CDF with freshly regenerated moment/transform shell prefix sums')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    run(**vars(args))


if __name__ == '__main__':
    main()
