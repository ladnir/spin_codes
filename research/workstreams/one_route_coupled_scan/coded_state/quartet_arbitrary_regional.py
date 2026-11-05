"""Compose arbitrary-pair edge bounds through complete restricted regions.

FLOAT only.1200s watchdog,<8GiB conservative arrays, single-thread Peach,
no files. All selected-count paths; not all counts or a whole-code cert.
"""
from collections import defaultdict
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_joint_regional as joint
import quartet_mixed_counts as mixed
from quartet_quota_potential import normalize,fit
from quartet_pair_regional import paired
from quartet_arbitrary_edges import coarse_cell,edge_caps


def tighten(matrix,cells,counts,labels,x,logv,growth,initial_v,bucket=64):
    start=time.monotonic();D=len(cells);groups=defaultdict(list)
    for source,cell in enumerate(cells+[(0,0,0,0,0,0)]):groups[coarse_cell(cell,bucket)].append(source)
    logv_all=np.r_[logv,initial_v];improved=0
    for index,(cell,sources) in enumerate(groups.items()):
        caps=edge_caps(cell,counts)
        for source in sources:
            normalized=caps+logv[None,:]-logv_all[source]-2*growth
            normalized+=x*(np.asarray(counts)[:,None]+labels[None,:]-816)
            with np.errstate(divide='ignore'):old=np.log2(matrix[source].reshape(len(counts),D))
            improved+=int(np.sum(normalized<old-1e-8))
            matrix[source]=mixed.positive_exp(np.minimum(old,normalized),'arbitrary_pair_edges').ravel()
        if index%8==0 or index+1==len(groups):
            print(dict(stage='arbitrary_pairs',coarse_group=index+1,groups=len(groups),
                       improved_edges=improved,seconds=time.monotonic()-start),flush=True)
    assert np.all(np.isfinite(matrix)) and max(matrix.sum(axis=1))<=1+1e-10
    return matrix


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(1200,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        counts=(0,312,384,392,407,408,424,512);joint.BIN=16;mixed.BIN=16;joint.states.cache_clear()
        base,scale,labels,offsets=joint.build(counts);cells=[];edge_ids=np.full((513,33,33),-1,dtype=int)
        for a in counts:
            for rl,rh,sl,sh in joint.states(a)[0]:
                edge_ids[a,rl//16,sl//16]=len(cells);cells.append((a,a,rl,rh,sl,sh))
        P,x,logv,growth,initial_v=normalize(base,scale,cells,edge_ids,labels,np.arange(len(cells)))
        B,pair_labels,targets=paired(P,cells,counts,offsets,labels,x,logv,growth,initial_v)
        del base,P
        B=tighten(B,cells,counts,labels,x,logv,growth,initial_v)
        result=fit(B,pair_labels,targets,cells,x,logv,growth,initial_v,iterations=4,steps=8,center=816)
        result.update(total_seconds=time.monotonic()-start,all_input_counts=False,
                      restricted_counts=list(counts),input_count_range=[0,512],two_update_blocks=True,
                      coupled_pairs='all pairs; arbitrary source buckets64, boundary buckets16',cells=len(cells))
        print(dict(arbitrary_paired_restricted_best=result),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
