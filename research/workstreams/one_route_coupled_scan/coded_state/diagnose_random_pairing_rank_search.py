"""Bounded route-and-pairing search for a positive structural rank gap.

Actual uniform regional routes and perfect matchings; no field matrices.
A positive gap is sufficient for an obstruction. A negative search is not
a certificate. Serial <=30s, <128MiB working arrays, no output files.
"""
import math
import shutil
import tempfile
import time
import numpy as np


def evaluate(route,partners,masks,block,rows=8192,parity=128,W=964689):
    states=np.repeat(masks,block,axis=0)
    previous=np.zeros_like(states);previous[1:]=states[:-1]
    live=previous|np.take_along_axis(previous,partners,axis=1)
    unions=previous|states;visible=live|states
    counts=np.bincount(route[visible],minlength=rows)
    J=int(np.count_nonzero(live&~unions));D=int(np.count_nonzero(previous&~states))
    area=int(np.count_nonzero(states));union=int(np.count_nonzero(unions))
    I=int(np.count_nonzero(previous&states));z=max(0,union-W)
    d=int(np.maximum(counts-parity,0).sum())
    gap=d-J-D-z
    assert gap==area-parity*rows+int(np.maximum(parity-counts,0).sum())-z
    assert int(np.count_nonzero(visible))==area+J+D
    return counts,dict(swaps=D,outside_crossings=J,shortened_dimension=d,
                       selected_output_zeros=z,rank_gap=gap,feasible_zero_count=z<=I,
                       structural_obstruction=gap>0 and z<=I,
                       parity_deficit=int(np.maximum(parity-counts,0).sum()),
                       rows_below_parity=int(np.count_nonzero(counts<parity)),
                       minimum_row_coverage=int(counts.min()))


def block_marginals(target,partners,masks,j,block):
    """Exact single-flip costs for coverage on target rows plus deaths.

    Death cost drops a linear cardinality term; it cancels in a swap.
    Returns the pair interaction needed to correct simultaneous swaps.
    """
    H=masks[j];start=j*block;end=start+block;L=len(target)
    if start:
        P=masks[j-1];live=P|P[partners[start]]
        cost=target[start].astype(np.int32)*(~live)-P.astype(np.int32)
    else:cost=target[0].astype(np.int32).copy()
    interior=partners[start+1:end]
    weights=target[start+1:end].astype(np.int32)
    pair_weights=weights+np.take_along_axis(weights,interior,axis=1)
    cost+=np.sum(pair_weights*(~H[interior]),axis=0,dtype=np.int32)
    last_partners=None;last_weights=None
    if end<L:
        N=masks[j+1];last_partners=partners[end]
        available=target[end].astype(np.int32)*(~N)
        last_weights=available+available[last_partners]
        cost+=last_weights*(~H[last_partners])-N.astype(np.int32)
    flip=(1-2*H.astype(np.int32))*cost
    return flip,interior,pair_weights,last_partners,last_weights


def improve_block(target,partners,masks,j,block):
    flip,ps,ws,lastp,lastw=block_marginals(target,partners,masks,j,block)
    H=masks[j];active=np.flatnonzero(H);inactive=np.flatnonzero(~H)
    aa=active[np.argsort(flip[active],kind='stable')[:8]]
    bb=inactive[np.argsort(flip[inactive],kind='stable')[:8]]
    best=(0,-1,-1)
    for a in aa:
        corrections=np.sum((ps[:,a,None]==bb[None,:])*ws[:,a,None],axis=0)
        if lastp is not None:corrections+=(lastp[a]==bb)*lastw[a]
        costs=flip[a]+flip[bb]+corrections
        k=int(np.argmin(costs));proposal=(int(costs[k]),int(a),int(bb[k]))
        if proposal<best:best=proposal
    if best[0]<0:
        H[best[1]]=False;H[best[2]]=True
        return True
    return False


def run():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--block',type=int,choices=(4,16,64))
    parser.add_argument('--independent-initial',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();rng=np.random.default_rng(20261027)
    rows=8192;R=256;n=512;L=16*R
    route=np.stack([rng.permutation(rows).astype(np.uint16) for _ in range(R)]).reshape(L,n)
    partners=np.empty((L,n),dtype=np.uint16)
    for t in range(L):
        order=rng.permutation(n);partners[t,order]=order[np.arange(n)^1]
    assert np.all(np.take_along_axis(partners,partners,axis=1)==np.arange(n))
    best=None;stopped=False
    for block in ((args.block,) if args.block else (64,16,4)):
        masks=np.zeros((L//block,n),dtype=bool);masks[:,:256]=True
        if args.independent_initial:
            for H in masks:
                H[:]=False;H[rng.permutation(n)[:256]]=True
        counts,result=evaluate(route,partners,masks,block)
        record=dict(result,block=block,sweep=-1)
        for sweep in range(12 if args.block else 6):
            selected=np.zeros(rows,dtype=bool)
            selected[np.argsort(counts,kind='stable')[:1974]]=True
            target=selected[route];accepted=0
            for j in range(len(masks)):
                for _ in range(8 if args.block else 4):
                    if not improve_block(target,partners,masks,j,block):break
                    accepted+=1
                if time.monotonic()-began>25:
                    stopped=True;break
            counts,result=evaluate(route,partners,masks,block)
            if result['rank_gap']>record['rank_gap']:
                record=dict(result,block=block,sweep=sweep,accepted_swaps=accepted)
            if stopped or result['structural_obstruction'] or not accepted:break
        print(dict(best=record,last=result,last_sweep=sweep,
                   lowest_1974_mean_coverage=float(np.sort(counts)[:1974].mean())),flush=True)
        if best is None or record['rank_gap']>best['rank_gap']:best=record
        if stopped or record['structural_obstruction']:break
    print(dict(seed=20261027,seconds=time.monotonic()-began,time_budget_reached=stopped,
               best_rank_gap=best['rank_gap'],whole_code_certificate=False,
               native_encoder_benchmark=False),flush=True)


if __name__=='__main__':run()
