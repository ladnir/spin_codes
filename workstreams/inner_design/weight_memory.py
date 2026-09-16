"""Pointwise weight-class envelope for lazy-refresh inners; binary64 Q1 only.

Unlike the three-state screen, live weight classes retain the expansion weight
through empty lazy epochs. A class coordinate is a mass envelope: coordinate
c at weight w dominates every state's mass by c/a_w. It does not assert that
the actual conditional state is uniform within a weight class.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import transvection as tv

HERE = Path(__file__).resolve().parent


def all_weights(columns, s):
    assert len(columns) <= 128
    halves = []
    for start in (0,64):
        generators = [sum(((columns[p]>>j)&1) << (p-start)
                          for p in range(start,min(start+64,len(columns)))) for j in range(s)]
        words = np.array([0],dtype=np.uint64)
        for row in generators:
            words = np.concatenate((words, words ^ np.uint64(row)))
        halves.append(np.bitwise_count(words))
    return halves[0]+halves[1]


def build_neighbors(columns, s):
    """Enumerate every state and single input exactly; retain Pareto pairs.

    Pair (minus,plus) for old weight v and destination weight w counts inputs
    with output weight v-1 or v+1. We take the maximum over destination states
    separately for each old class; this is conservative, not a lumpability claim.
    """
    weights = all_weights(columns,s)
    spectrum = {int(w):int(n) for w,n in zip(*np.unique(weights,return_counts=True)) if w}
    assert sum(spectrum.values()) == (1<<s)-1
    levels = sorted(spectrum)
    indices = np.zeros(len(columns)+1,dtype=np.uint8)
    for i,w in enumerate(levels): indices[w]=i
    states = np.arange(1<<s,dtype=np.uint32)
    counts = np.zeros((len(states),len(levels),2),dtype=np.uint8)
    for p,col in enumerate(columns):
        old = states ^ col
        bit = np.bitwise_count(old & col) & 1
        keep = old != 0
        # Index 0 is the minus (bit=1) count; index 1 the plus count.
        counts[states[keep], indices[weights[old[keep]]], 1-bit[keep]] += 1
    pairs = []
    for old_index,v in enumerate(levels):
        destinations=[]
        for w in levels:
            possibilities = np.unique(counts[weights==w,old_index,:],axis=0)
            frontier = [p.tolist() for p in possibilities
                        if not any(np.all(other>=p) and np.any(other>p) for other in possibilities)]
            destinations.append(frontier)
        pairs.append(destinations)
    cw = tv.cancellation_weights(columns)
    cancellation = [[int(weights[col]), cw[p]] for p,col in enumerate(columns)]
    return dict(s=s,t=len(columns),spectrum=spectrum,levels=levels,
                lazy_neighbors=pairs,cancellation=cancellation,
                checked_state_input_pairs=len(states)*len(columns))


def transfers(record,lambdas,rounds,active_mode='forget'):
    lam=np.asarray(lambdas)
    s,t=record['s'],record['t'];den=(1<<s)-1
    spectrum={int(w):n for w,n in record['spectrum'].items()}
    levels=record['levels'];n=len(levels)+2
    eps=0. if rounds is None else 2.**-rounds
    le=-np.inf if eps==0 else math.log(eps);lf=math.log1p(-eps)
    f0=-lam[:,None]*np.array(levels)
    f1=np.empty_like(f0)
    for i,w in enumerate(levels):
        f1[:,i]=-lam*(t-1) if w==t else np.logaddexp(
            math.log(w/t)-lam*(w-1),math.log1p(-w/t)-lam*(w+1))
    d0=np.max(f0,axis=1);d1=np.max(f1,axis=1)
    zero=np.full((len(lam),n,n),-np.inf);one=np.full_like(zero,-np.inf)
    zero[:,0,0]=0;one[:,0,1]=-lam
    for matrix,f,d in ((zero,f0,d0),(one,f1,d1)):
        matrix[:,1,1]=le+d
        for j,w in enumerate(levels):
            matrix[:,1,j+2]=lf+d+math.log(spectrum[w]/den)
            for i,v in enumerate(levels):
                matrix[:,i+2,j+2]=lf+f[:,i]+math.log(spectrum[w]/den)
    for i,v in enumerate(levels):
        zero[:,i+2,i+2]=np.logaddexp(zero[:,i+2,i+2],le+f0[:,i])
        lazy_cancel=[y for w,y in record['cancellation'] if w==v]
        if lazy_cancel:
            mass=np.logaddexp.reduce(-lam[:,None]*np.array(lazy_cancel),axis=1)-math.log(t*spectrum[v])
        else: mass=np.full(len(lam),-np.inf)
        one[:,i+2,0]=np.logaddexp(le+mass,lf+f1[:,i]-math.log(den))
        if active_mode=='forget':
            # Only a nonempty lazy epoch forgets the class. Its total emitted
            # mass is bounded exactly by f1; do not inflate all destination caps.
            one[:,i+2,1]=le+f1[:,i]
            continue
        assert active_mode=='neighbors'
        for j,w in enumerate(levels):
            largest=np.full(len(lam),-np.inf)
            for minus,plus in record['lazy_neighbors'][i][j]:
                left=-np.inf if minus==0 else math.log(minus)-lam*(v-1)
                right=-np.inf if plus==0 else math.log(plus)-lam*(v+1)
                largest=np.maximum(largest,np.logaddexp(left,right))
            value=le+largest+math.log(spectrum[w]/(t*spectrum[v]))
            one[:,i+2,j+2]=np.logaddexp(one[:,i+2,j+2],value)
    cd=-lam*min(y for w,y in record['cancellation'])-math.log(t)
    one[:,1,0]=np.logaddexp(le+cd,lf+d1-math.log(den))
    return zero,one


def product(a,b):
    result=a[:,:,0,None]+b[:,None,0,:]
    for i in range(1,a.shape[1]):
        np.logaddexp(result,a[:,:,i,None]+b[:,None,i,:],out=result)
    return result


def regions(zero,one,epochs):
    rz=np.full_like(zero,-np.inf)
    for i in range(zero.shape[1]):rz[:,i,i]=0
    ra=np.full_like(one,-np.inf)
    while epochs:
        if epochs&1:
            ra=np.logaddexp(product(ra,zero),product(rz,one));rz=product(rz,zero)
        epochs>>=1
        if epochs:
            one=np.logaddexp(product(one,zero),product(zero,one));zero=product(zero,zero)
    return rz,ra


def coefficients(rz,ra,length):
    current=np.full((len(rz),length+1,rz.shape[1]),-np.inf);current[:,0,0]=0
    def vm(v,m):
        out=v[:,:,0,None]+m[:,None,0,:]
        for i in range(1,v.shape[2]):np.logaddexp(out,v[:,:,i,None]+m[:,None,i,:],out=out)
        return out
    for n in range(length):
        updated=np.full_like(current,-np.inf)
        updated[:,:n+1]=vm(current[:,:n+1],rz)
        updated[:,1:n+2]=np.logaddexp(updated[:,1:n+2],vm(current[:,:n+1],ra))
        current=updated
    return np.logaddexp.reduce(current,axis=2)-np.array([math.log(math.comb(length,w)) for w in range(length+1)])


def main():
    # Recompute the exact neighbor envelope; no unverified cached table input.
    spectrum,_,sources=tv.fixed.load_inner()
    path=tv.ROOT/'workstreams/bare_bch_rm2sub/generated/manifest.json'
    columns=json.loads(path.read_text())['t128_s19']['columns']
    record=build_neighbors(columns,19)
    assert record['spectrum']==spectrum
    neighbors=HERE/'WEIGHT_NEIGHBORS.json'
    neighbors.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('Exact single-input neighbors enumerated.',flush=True)
    outer={w:n for w,n in enumerate(tv.smaller_outer.spectrum()) if w and n}
    tilts=np.arange(-1020,-879,dtype=float)/100
    rows=32768;results=[]
    for active_mode,rounds in [('forget',r) for r in (None,1,2,3,4)]+[('neighbors',1)]:
        rz,ra=regions(*transfers(record,np.exp(tilts),rounds,active_mode),rows//128)
        moments=coefficients(rz,ra-math.log(rows//128),128)
        for delta in (tv.Fraction(33,200),tv.Fraction(19,100)):
            cutoff=128*rows*delta.numerator//delta.denominator
            values=np.minimum(0.,moments+cutoff*np.exp(tilts)[:,None])
            ix=np.argmin(values,axis=0);best=values[ix,np.arange(129)]
            terms=[math.log(rows*n)+best[w] for w,n in outer.items()]
            dominant=list(outer)[int(np.argmax(terms))]
            row=dict(active_mode=active_mode,rounds=rounds,distance_target=str(delta),
                     q1_margin_bits=-float(np.logaddexp.reduce(terms))/math.log(2),
                     dominant_outer_weight=dominant,dominant_log_surprisal=float(tilts[ix[dominant]]),
                     witness_at_grid_edge=bool(ix[dominant] in (0,len(tilts)-1)))
            results.append(row);print(json.dumps(row),flush=True)
    sources += [Path(__file__),Path(tv.__file__),Path(tv.smaller_outer.__file__),path,neighbors]
    payload=dict(status='BINARY64_Q1_ONLY_NOT_A_CERTIFICATE',outer=[128,32,32],
                 message_bits=1<<20,t=128,s=19,q1_log_surprisals=tilts.tolist(),
                 source_sha256={p.relative_to(tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources},results=results)
    (HERE/'WEIGHT_MEMORY_Q1.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
