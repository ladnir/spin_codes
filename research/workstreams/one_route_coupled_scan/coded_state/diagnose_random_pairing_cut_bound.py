"""Bounded fixed-target cut diagnostic on sampled full-size geometry.

Verified integer bounds for processed windows only. No probability bound
over setup or all target sets, and no distance certificate. No files.
"""
import argparse
import shutil
import tempfile
import time
import numpy as np
from random_pairing_cut_bound import cut_bound
from diagnose_random_pairing_rank_search import evaluate


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true')
    parser.add_argument('--window',type=int,choices=(16,32,64),default=32)
    parser.add_argument('--symmetric',action='store_true')
    parser.add_argument('--scale',type=int,choices=(16,64,256),default=16)
    parser.add_argument('--mu',type=int)
    parser.add_argument('--start-window',type=int,default=0)
    args=parser.parse_args()
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();rng=np.random.default_rng(20261027)
    M=8192;L=4096;n=512
    route=np.stack([rng.permutation(M).astype(np.uint16) for _ in range(256)]).reshape(L,n)
    partners=np.empty((L,n),dtype=np.uint16)
    for t in range(L):
        order=rng.permutation(n);partners[t,order]=order[np.arange(n)^1]
    masks=np.zeros((1,n),dtype=bool);masks[0,:256]=True
    counts,_=evaluate(route,partners,masks,L)
    selected=np.zeros(M,dtype=bool);selected[np.argsort(counts,kind='stable')[:1974]]=True
    target=selected[route]
    assert 0<=args.start_window<L//args.window
    starts=list(range(args.start_window*args.window,L,args.window)) if args.full else [0,1024,2048,3072]
    total=0;completed=0;stopped=False;stop_reason=None
    for start in starts:
        best=None;trials=[]
        for mu in ((args.mu,) if args.mu is not None else (3,4,5)):
            if shutil.disk_usage(tempfile.gettempdir()).free<20*(1<<30):
                stopped=True;stop_reason='disk guard';break
            try:
                result=cut_bound(target[start:start+args.window],partners[start:start+args.window],mu,
                                 scale=args.scale,symmetric=args.symmetric)
            except MemoryError:
                stopped=True;stop_reason='allocation failure';break
            trials.append(result)
            if best is None or result['lower_scaled']>best['lower_scaled']:best=result
        if best is None:break
        total+=max(0,best['lower_scaled']);completed+=1
        if not args.full:print(dict(start=start,window=args.window,trials=trials),flush=True)
        if time.monotonic()-began>25:
            stopped=True;stop_reason='time budget'
        if stopped:break
    # Unprocessed windows have the valid lower bound zero, not an extrapolation.
    whole_partition=args.full and args.start_window==0 and completed==len(starts)
    obstruction_threshold=128*1974-83887
    print(dict(seed=20261027,window=args.window,processed_windows=completed,
               first_window=starts[0]//args.window,
               last_window=starts[completed-1]//args.window if completed else None,
               processed_lower_scaled=total,scale=args.scale,
               symmetric_deaths=args.symmetric,
               processed_cost_lower=total/args.scale,whole_partition_processed=whole_partition,
               fixed_target_rank_gap_upper=obstruction_threshold-total/args.scale,
               time_budget_reached=stop_reason=='time budget',stop_reason=stop_reason,
               seconds=time.monotonic()-began,
               all_target_sets_covered=False,whole_code_certificate=False),flush=True)


if __name__=='__main__':main()
