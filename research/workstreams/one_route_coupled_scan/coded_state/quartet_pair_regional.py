"""Complete restricted count paths, improving every aligned full-input pair.

FLOAT diagnostic only. Serial Peach,1200s,<8GiB conservative owned arrays,
no generated files. One pair contains two actual updates, not a new code.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_joint_regional as joint
import quartet_mixed_counts as mixed
from quartet_full_pair import prepare_first,prepare_second,full_first,full_second
from quartet_quota_potential import normalize,fit


def paired(P,cells,counts,offsets,labels,x,logv,growth,initial_v):
    start=time.monotonic();D=len(cells);E=len(counts)*D
    assert (4*(D+1)*E+2*8193*(D+1+E))*8<8*(1<<30)
    matrix=np.empty((D+1,E));where={a:i for i,a in enumerate(counts)}
    pair_labels=np.concatenate([labels+a for a in counts]);targets=np.tile(np.arange(D),len(counts))
    for i,a in enumerate(counts):
        mid=slice(offsets[i],offsets[i+1])
        matrix[:,i*D:(i+1)*D]=np.maximum(P[:,mid]@P[mid,:D],2**mixed.FLOOR)
    logv_all=np.r_[logv,initial_v];improved=0
    def intersect(source,first,last,raw):
        nonlocal improved
        i=where[first];j=where[last];columns=np.arange(offsets[j],offsets[j+1])
        cap=np.full(len(columns),-np.inf);lookup=joint.states(last)[1]
        for key,value in raw.items():
            assert key[0]<lookup.shape[0] and key[1]<lookup.shape[1] and lookup[key]>=0
            cap[lookup[key]]=value
        cap+=logv[columns]-logv_all[source]-2*growth+x*(first+last-816)
        edge=i*D+columns
        with np.errstate(divide='ignore'):old=np.log2(matrix[source,edge])
        improved+=int(np.sum(cap<old-1e-8))
        matrix[source,edge]=mixed.positive_exp(np.minimum(old,cap),'paired_edge')
    for source,cell in enumerate(cells+[(0,0,0,0,0,0)]):
        first=prepare_first(cell);second=prepare_second(cell)
        for l in counts:
            intersect(source,512,l,full_first(cell,l,prepared=first))
            intersect(source,l,512,full_second(cell,l,prepared=second))
        if source%128==0 or source==D:
            print(dict(stage='coupled_full_pairs',source=source,states=D,improved_edges=improved,
                       seconds=time.monotonic()-start),flush=True)
    assert np.all(np.isfinite(matrix)) and np.max(matrix.sum(axis=1))<=1+1e-10
    return matrix,pair_labels,targets


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(1200,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        counts=(0,312,384,392,407,408,424,512);joint.BIN=16;mixed.BIN=16;joint.states.cache_clear()
        base,scale,labels,offsets=joint.build(counts);cells=[];edge_ids=np.full((513,33,33),-1,dtype=int)
        for i,a in enumerate(counts):
            for j,(rl,rh,sl,sh) in enumerate(joint.states(a)[0]):
                edge_ids[a,rl//16,sl//16]=len(cells);cells.append((a,a,rl,rh,sl,sh))
        targets=np.arange(len(cells))
        P,x,logv,growth,initial_v=normalize(base,scale,cells,edge_ids,labels,targets)
        B,paired_labels,paired_targets=paired(P,cells,counts,offsets,labels,x,logv,growth,initial_v)
        del base,P
        result=fit(B,paired_labels,paired_targets,cells,x,logv,growth,initial_v,iterations=4,steps=8,center=816)
        result.update(total_seconds=time.monotonic()-start,all_input_counts=False,
                      restricted_counts=list(counts),input_count_range=[0,512],two_update_blocks=True,
                      coupled_pairs='one or both counts512',cells=len(cells))
        print(dict(paired_restricted_best=result),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
