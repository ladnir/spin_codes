"""Regenerate selected fixed-occupancy support covers for shared-row GF16.

Only the explicitly listed occupancies are covered. This is not a full-code
certificate; in particular, the independent-row dense cover does not apply.
"""
import argparse
import json
import hashlib
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from flint import arb, ctx
import sparse_cover
import occupancy_birth_classes
from shared_support import shared_counts


def run(updates, occupancies, threshold, tilts, precision, max_splits, target_bits, output=None, probe_supports=(), refined_counts=False, proposal_buffer_bits=2, joint_return_through=None, lazy_density_through=None, coupled_counts=False, joint_counts=False, wide_counts=False, count_witnesses=()):
    if (not updates or any(type(r) is not int or r not in (2,3,4) for r in updates)
            or not occupancies or len(set(occupancies)) != len(occupancies)
            or any(type(q) is not int or not 1 <= q <= 2048 for q in occupancies)
            or precision < 128 or max_splits < 0
            or type(target_bits) is not int or target_bits < 1
            or type(threshold) is not int or not 0 <= threshold < 1 << 21
            or not tilts or any(Q(t) <= 0 for t in tilts)
            or not 0 <= proposal_buffer_bits <= 64
            or any(type(u) is not int or not 38 <= u <= 256 for u in probe_supports)):
        raise ValueError('valid update counts, distinct occupancies, and search limits required')
    for value, limit in ((joint_return_through, 4), (lazy_density_through, 32)):
        if value is not None and (type(value) is not int or not 0 <= value <= limit):
            raise ValueError('valid integer conditional GF refinement cutoffs required')
    ctx.prec = precision
    counts, _ = shared_counts(refined_counts, coupled_counts, joint_counts, wide_counts, count_witnesses)
    results = []
    for r in updates:
        args = sparse_cover.build_args(max(occupancies),tilts,precision,max_splits,target_bits,None,r)
        args.exact_feedback = True
        args.analytic_gradient = True
        args.probe_supports = list(probe_supports)
        args.proposal_buffer_bits = proposal_buffer_bits
        args.joint_return_through = joint_return_through
        args.lazy_density_through = lazy_density_through
        operators = occupancy_birth_classes.build_operators(args)
        terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
        for q in occupancies:
            args.groups = q
            print('SHARED GF16', 'PROBE' if probe_supports else 'COVER',
                  'updates',r,'occupancy',q,'cutoff',threshold,flush=True)
            upper = sparse_cover.sparse.cover(args,operators,{'1':counts},terminal,cutoff=threshold)
            results.append(dict(updates=r,occupancy=q,upper=None if upper is None
                                else [int(v) for v in upper.upper().man_exp()]))
            if output:
                record = dict(schema='shared-gf16-sparse-screen-1',precision=precision,
                              updates=updates,occupancies=occupancies,threshold=threshold,
                              tilts=tilts,max_splits=max_splits,target_bits=target_bits,results=results,
                              probe_supports=list(probe_supports),
                              refined_counts=refined_counts,
                              coupled_counts=coupled_counts,
                              joint_counts=joint_counts,
                              wide_counts=wide_counts,
                              count_witnesses=[dict(path=str(p), sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())
                                               for p in count_witnesses],
                              joint_return_through=joint_return_through,
                              lazy_density_through=lazy_density_through,
                              proposal_buffer_bits=proposal_buffer_bits,
                              note='Fresh shared-row group CDF and GF operators. Listed occupancies only; not a full-code certificate.')
                output.parent.mkdir(parents=True,exist_ok=True)
                output.write_text(json.dumps(record,indent=2)+'\n')
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--updates',type=int,nargs='+',default=[2,4])
    parser.add_argument('--occupancies',type=int,nargs='+',default=[2,4,8,16])
    parser.add_argument('--threshold',type=int,default=209715)
    parser.add_argument('--tilts',nargs='+',default=['.00032','.001','.0032','.008','.016','.032','.064','.096'])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--max-splits',type=int,default=64)
    parser.add_argument('--target-bits',type=int,default=48,
                        help='Positive per-occupancy proof budget; not the aggregate code margin')
    parser.add_argument('--probe-supports',type=int,nargs='+',default=[],help='Floating selected-support diagnostics only, not a support cover')
    parser.add_argument('--refined-counts',action='store_true',help='Regenerate positive shortening and dual-moment support bounds')
    parser.add_argument('--coupled-counts',action='store_true',help='Also alternate containment and dual bounds; requires --refined-counts')
    parser.add_argument('--joint-counts',action='store_true',help='Add exact-verified joint CDF constraints; requires --coupled-counts')
    parser.add_argument('--wide-counts',action='store_true',help='Also check wider CDF constraints through support 192; requires --joint-counts')
    parser.add_argument('--count-witnesses',type=Path,nargs='+',default=[],help='Replay saved exact dual multipliers against fresh counting constraints; requires --coupled-counts')
    parser.add_argument('--joint-return-through',type=int,choices=(0,1,2,3,4),help='Retain weighted returns through this local packet occupancy')
    parser.add_argument('--lazy-density-through',type=int,help='Refine conditional lazy-state densities through this local occupancy (0..32)')
    parser.add_argument('--proposal-buffer-bits',type=float,default=2,help='Search stopping buffer only; the directed result must still meet target-bits')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    run(args.updates,args.occupancies,args.threshold,args.tilts,args.precision,args.max_splits,args.target_bits,args.output,args.probe_supports,args.refined_counts,args.proposal_buffer_bits,args.joint_return_through,args.lazy_density_through,args.coupled_counts,args.joint_counts,args.wide_counts,args.count_witnesses)
