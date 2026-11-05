"""Uniform marked arbitrary-pair edge bounds; local diagnostic, not a cert."""
import math
import time
import numpy as np
import quartet_joint_regional as joint
import quartet_mixed_counts as mixed
from quartet_arbitrary_pair import moment_bounds
from quartet_full_pair import logsum


MARKERS=tuple((r,s) for r in (-.5,0.,.5) for s in (-.5,0.,.5))


def coarse_cell(cell,bucket=64):
    a,ahi,rl,rh,sl,sh=cell;assert a==ahi
    r=rl//bucket*bucket;s=sl//bucket*bucket
    assert rh<=min(a,r+bucket-1) and sh<=min(a-r,s+bucket-1)
    return a,a,r,min(a,r+bucket-1),s,min(a-r,s+bucket-1)


def edge_caps(cell,counts,markers=MARKERS):
    boxes=[joint.states(l)[0] for l in counts]
    offsets=np.cumsum([0]+[len(b) for b in boxes]);D=int(offsets[-1])
    caps=np.full((len(counts),D),np.inf)
    for r,s in markers:
        bounds=moment_bounds(cell,counts,counts,R=2**r,S=2**s)
        for j,l in enumerate(counts):
            endpoints=np.array([-r*(rl if r>=0 else rh)-s*(sl if s>=0 else sh)
                                for rl,rh,sl,sh in boxes[j]])
            for i,first in enumerate(counts):
                caps[i,offsets[j]:offsets[j+1]]=np.minimum(caps[i,offsets[j]:offsets[j+1]],bounds[first,l]+endpoints)
    return caps


def diagnostic():
    start=time.monotonic();joint.BIN=16;joint.states.cache_clear();counts=(312,408,512)
    base,scale,labels,offsets=joint.build(counts);where={a:i for i,a in enumerate(counts)}
    for a,rl,sl in ((312,64,16),(408,192,96),(512,0,0)):
        i=where[a];boxes,lookup=joint.states(a);j=lookup[rl//16,sl//16]
        cell=(a,a,*boxes[j]);coarse=coarse_cell(cell);source=offsets[i]+j
        caps=edge_caps(coarse,counts)
        for first,last in ((312,408),(408,312),(408,408),(512,312)):
            mi=where[first];oi=where[last];ms=slice(offsets[mi],offsets[mi+1]);dest=slice(offsets[oi],offsets[oi+1])
            left=base[source,ms];right=base[ms,dest];ls=max(left);rs=np.max(right)
            row=(left/ls)@(right/rs)
            # Preserve an upper representation even for positive underflow.
            positive=(left>0).astype(int)@(right>0).astype(int)>0
            row=np.where(positive,np.maximum(row,2**mixed.FLOOR),0.)
            with np.errstate(divide='ignore'):old=2*scale+math.log2(ls)+math.log2(rs)+np.log2(row)
            new=np.minimum(old,caps[mi,dest])
            print(dict(source=cell,cap_cell=coarse,counts=[first,last],ordinary=logsum(old),
                       intersected=logsum(new),improved_edges=int(np.sum(new<old-1e-8))),flush=True)
    print(dict(seconds=time.monotonic()-start,local_only=True),flush=True)
