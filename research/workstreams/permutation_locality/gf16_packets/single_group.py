"""Exact support-placement transfer for one active four-row group.

This is an occupancy-one bound, not a full-code certificate. The local
state envelopes remain upper bounds; only their placement average is exact.
"""
import argparse
import json
from pathlib import Path
from fractions import Fraction as Q
from flint import arb, arb_mat, ctx
import sparse_cover
import occupancy_rank

up = occupancy_rank.up


def support_moments(inactive, active, regions=256, *, matrix=arb_mat,
                    scalar=arb, rounding=up):
    """Average ordered products over all supports of each cardinality."""
    size=inactive.nrows()
    if (type(regions) is not int or regions<0 or size<1
            or inactive.ncols()!=size or active.nrows()!=size or active.ncols()!=size):
        raise ValueError('equal square matrices and nonnegative region count required')
    values=[matrix([[1]+[0]*(size-1)])]
    for r in range(1,regions+1):
        following=[]
        for j in range(r+1):
            row=matrix(1,size)
            if j<r:row+=values[j]*inactive*(scalar(r-j)/scalar(r))
            if j:row+=values[j-1]*active*(scalar(j)/scalar(r))
            following.append(matrix([[rounding(row[0,k]) for k in range(size)]]))
        values=following
    return [rounding(sum((row[0,k] for k in range(size)),scalar(0))) for row in values]


def fold_cdf(counts, weights, *, scalar=arb, rounding=up):
    """Use a nonincreasing majorant, not CDF differences as shell caps."""
    if (len(counts)!=len(weights) or not counts or counts[0]!=0
            or any(c<0 for c in counts) or any(a>b for a,b in zip(counts,counts[1:]))
            or any(w<0 for w in weights)):
        raise ValueError('nonzero-message CDF and nonnegative weights required')
    tail=[];largest=scalar(0)
    for value in reversed(weights):
        largest=max(largest,value);tail.append(largest)
    tail.reverse()
    return rounding(sum(((counts[u]-counts[u-1])*tail[u] for u in range(1,len(counts))),scalar(0)))


def run(threshold, tilts, precision=256, output=None,updates=2):
    if not 0<=threshold<(1<<21) or not tilts or any(Q(t)<=0 for t in tilts) or precision<128:
        raise ValueError('valid threshold, positive tilts and precision >=128 required')
    if type(updates) is not int or updates not in (2,3,4):raise ValueError('two through four inner updates required')
    ctx.prec=precision
    sparse=sparse_cover.sparse
    counts=sparse.integer_cdf(sparse.weighted_cdf_upper(sparse.authenticated_caps(),1<<128,full_weight=Q(1)))
    args=sparse_cover.build_args(1,tilts,precision,0,40,None,updates)
    operators=occupancy_rank.build_operators(args)
    best=[arb(1)]*257;choices=[None]*257
    for tilt in tilts:
        exact=operators[tilt,'1'][0]
        weights=support_moments(exact[0],exact[1])
        multiplier=(occupancy_rank.aq(Q(tilt))*threshold).exp()
        for u,value in enumerate(weights):
            value=up(value*multiplier)
            if value<best[u]:best[u]=value;choices[u]=tilt
        upper=up(2048*fold_cdf(counts,best))
        print('ONE GROUP exact placement through tilt',tilt,'log2 upper',upper.log()/arb(2).log(),flush=True)
    print('ONE GROUP complete support range, cutoff',threshold,'margin',-upper.log()/arb(2).log(),flush=True)
    print('Other occupancies are not covered.',flush=True)
    if output:
        record=dict(schema='gf16-single-group-run-1',threshold=threshold,precision=precision,updates=updates,
                    tilts=tilts,upper=[int(x) for x in upper.upper().man_exp()],
                    support_choices=choices,
                    support_probability_uppers=[[int(x) for x in v.upper().man_exp()] for v in best],
                    note='Occupancy one only. Regenerate with this program; not a full-code certificate.')
        output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(record,indent=2)+'\n')
    return upper


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threshold',type=int,default=209715)
    parser.add_argument('--tilts',nargs='+',default=['.00024','.00028','.00032','.0004','.0032','.016','.032','.064','.096'])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--updates',type=int,choices=(2,3,4),default=2)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    run(args.threshold,args.tilts,args.precision,args.output,args.updates)


if __name__=='__main__':main()
