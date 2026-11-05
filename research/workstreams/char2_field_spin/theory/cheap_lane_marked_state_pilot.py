"""FLOAT complementary-branch pilot with actual alternating source profiles.

Only consecutive sets of256 source groups are considered:4096 sets have
regional histograms(128 full,0 single),(127 full,2 single), and4096 have
the reverse order. No artificial repeated marginal is used.

Raw per-fixed-destination scalar caps are averaged at alpha1, with actual
support filters, deterministic gamma minima, and the original joint-live
cap. Raw matrices include destination multiplicity D_j=C(16,j). Pair-only
packet occupancies additionally use the frozen shared-point singleton
extinction cap. Asymmetric singleton occupancies retain the generic cap.

Insert a quiet mark t on00, or a birth mark b on0j,j>0, BEFORE pooling
internal masks and taking alpha. For coefficient markers u,w, let Mbar
and Tbar be the raw packet matrix and output-only vector averaged with
weights coefficient*u**f*w**s/U**32, U=1+u+w. Packet singleton degrees>2
are set to zero: they cannot occur in either target regional coefficient.
The full route mass U**32 is retained in Jensen, so this deletion is safe.
All full-pair tails and all singleton degrees0,1,2 are included.

A tile of p=16 packets has normalized Jensen transfer
 D_j**(1-alpha)*(Mbar**p)_hj**alpha.
Its terminal vector is (Mbar**(p-1)*Tbar)_h**alpha, with no final mask,
birth mark, quiet mark, or feedback. Powers here are entrywise only where
alpha is written. Matrix powers sum all internal masks before alpha.
The regional Cauchy factor is exp(H), H=4096logU-f logu-s logw-log M4096.
It is NOT raised to alpha. Eight tiles form each128-packet region.
After factoring exp(128*(H0+H1)), the two orientation chains are
 e0*(R_i R_(1-i))**127*R_i*V_(1-i), i=0,1.
Their positive sum retains exactly32767 regular transitions.

Class II uses b>=1 and B>=513, with correction b**(-alpha*513).
Class III uses0<t<=1 and Q<=24575, with correction t**(-alpha*24575).
The latter drops the birth-count restriction, a safe enlargement. Both
also impose output weight<=floor(.35*2**21). Source and output factors
are (eta*beta**256*z**(-L))**alpha, eta=256**2/255**2,
beta=256**160/255**32. The4096 multiplicity per orientation is outside.
An additional raw row mark a**ceil(h/4), terminal included, enforces the
necessary global floor S=sum ceil(h_t/4)<=61503. Indeed W>=13S-256r;
the correction is a**(-alpha*61503), with0<a<=1. It uses actual masks,
not the auxiliary scalar state. a=1 recovers the unmarked relaxation.

Four fixed(z,alpha) anchors, two branches, single-start marker optimization.
Finite positive log-matrix products, never Perron surrogates. This is a
selected-source-class FLOAT diagnostic, not a whole-code certificate.
One48-second invocation,20GiB free, under128MiB estimated peak. No data,
cache, frozen-source edits, broad scans, or additional tile runs.
"""

from collections import defaultdict
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
PIN="4526522c47c9920482f8990bb104243eada12318774ee2d0dbb68fda3704ab36"
LN2=math.log(2)
LIMIT=35*(1<<21)//100
DEADLINE=0.
np=None
logsumexp=None


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("48-second deadline exhausted; no automatic retry")


def imports():
    path=Path(__file__).with_name("cheap_lane_grouped_masks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen grouped-mask helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_marked_grouped",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    singleton=module.load_singleton()
    paired=singleton.load_block().load_paired(DEADLINE)
    filtered,_=paired.imports()
    return singleton,filtered


def geometry():
    """Positive exact dyadic counts; retain singleton total<=2 only."""
    half=[]
    for f in range(17):
        for s in range(min(2,16-f)+1):
            for left in range(s+1):
                coefficient=Fraction(math.comb(16,f)*math.comb(16-f,s)*math.comb(s,left),2**s)
                half.append((f,s,(f+left,f+s-left),coefficient))
    packet=defaultdict(Fraction)
    for f,s,shape,w in half:
        for f2,s2,shape2,w2 in half:
            if s+s2<=2:
                packet[(tuple(sorted(shape+shape2)),f+f2,s+s2)]+=w*w2
    shapes=sorted({key[0] for key in packet}); index={shape:i for i,shape in enumerate(shapes)}
    records=sorted(packet)
    checked=0
    for s in range(3):
        for f in range(33-s):
            got=sum(value for (shape,full,single),value in packet.items() if full==f and single==s)
            assert got==math.comb(32,f)*math.comb(32-f,s)
            checked+=1
    return (np.asarray(shapes,dtype=np.int16),
            np.asarray([index[key[0]] for key in records],dtype=np.int32),
            np.asarray([(key[1],key[2]) for key in records],dtype=float),
            np.asarray([float(packet[key]) for key in records]),records,checked)


def raw_tables(shapes,index,records,z,singleton,filtered):
    matrix=terminal=None
    for gamma in (0.,.5,1.):
        guard()
        local,final=filtered.local_tables(shapes,z,1.,gamma,1.)
        for j in range(1,17):
            factor=min(1.,256/255*255.**(-(16-j)))
            local[:,:,:,j]=np.minimum(local[:,:,:,j],final[:,:,:,0]*factor*math.comb(16,j))
        m=local[:,0]+local[:,2]; t=final[:,0]+final[:,2]
        matrix=m if matrix is None else np.minimum(matrix,m)
        terminal=t if terminal is None else np.minimum(terminal,t)
    matrix=matrix[index].copy(); terminal=terminal[index].copy()
    halves=singleton.half_tables(z,1.)
    refined=0
    for i,(shape,_,s) in enumerate(records):
        if s==0:
            assert shape[0]==shape[1] and shape[2]==shape[3]
            for h in range(1,17):
                matrix[i,h,0]=min(matrix[i,h,0],singleton.packet_cap(halves,shape[0],shape[2],h))
            refined+=1
    assert np.all(matrix>=0) and np.isfinite(matrix).all()
    assert np.all(terminal>=0) and np.isfinite(terminal).all()
    return matrix,terminal,refined


def tolog(value):
    out=np.full(value.shape,-math.inf)
    np.log(value,out=out,where=value>0)
    return out


def mm(a,b):
    return logsumexp(a[:,:,None]+b[None,:,:],axis=1)


def power(a,n):
    out=np.full(a.shape,-math.inf); np.fill_diagonal(out,0.)
    while n:
        if n&1:
            out=mm(out,a)
        n//=2
        if n:
            a=mm(a,a)
    return out


def double_pair(a,b,times):
    for _ in range(times):
        b,a=mm(a,b),mm(a,a)
    return a,b


def region(logm,logt,alpha,logd):
    a,b=double_pair(logm,logt,4)
    a=alpha*a+(1-alpha)*logd[None,:]
    b=alpha*b
    return double_pair(a,b,3)


def chain(regions,tails):
    values=[]
    for i in (0,1):
        total=mm(mm(power(mm(regions[i],regions[1-i]),127),regions[i]),tails[1-i])
        values.append(float(total[0,0]))
    return float(logsumexp(values)),values


def tests():
    """Explicit mask and route enumeration with nonuniform endpoint counts."""
    d=(1,2,1); categories=(0,1,1,2); alpha=.3
    caps=[]; ends=[]
    mark=Fraction(3,4)
    for occupied in (0,1):
        k=[[Fraction(2+i+j+occupied,11+2*i+3*j+occupied) for j in range(3)] for i in range(3)]
        k[0][0]*=mark
        caps.append(k)
        ends.append([[Fraction(1+i+occupied,9+2*i+occupied)] for i in range(3)])
    odds=Fraction(2,3); U=1+odds
    m=[[sum(odds**r*caps[r][i][j] for r in (0,1))*d[j]/U for j in range(3)] for i in range(3)]
    t=[[sum(odds**r*ends[r][i][0] for r in (0,1))/U] for i in range(3)]
    direct=[[Fraction(0) for _ in range(3)] for _ in range(3)]
    terminal=[[Fraction(0)] for _ in range(3)]
    powered=np.zeros((3,3)); powered_terminal=np.zeros((3,1))
    for a,b in itertools.product((0,1),repeat=2):
        route=odds**(a+b)/U**2
        for i in range(3):
            for j in range(3):
                raw=sum(caps[a][i][mid]*caps[b][mid][j] for mid in categories)
                direct[i][j]+=route*raw*d[j]
                powered[i,j]+=float(route)*float(raw)**alpha*d[j]
            raw=sum(caps[a][i][mid]*ends[b][mid][0] for mid in categories)
            terminal[i][0]+=route*raw
            powered_terminal[i,0]+=float(route)*float(raw)**alpha
    ma=np.asarray(m,dtype=float); ta=np.asarray(t,dtype=float)
    assert np.allclose(ma@ma,np.asarray(direct,dtype=float),rtol=1e-14,atol=0)
    assert np.allclose(ma@ta,np.asarray(terminal,dtype=float),rtol=1e-14,atol=0)
    bound=(ma@ma)**alpha*np.asarray(d,dtype=float)[None,:]**(1-alpha)
    assert np.all(powered<=bound*(1+1e-14))
    assert np.all(powered_terminal<=(ma@ta)**alpha*(1+1e-14))
    # Honest terminal: two packets have exactly one quiet/birth mark.
    bare=[[Fraction(1,2),Fraction(1,3)],[Fraction(1,4),Fraction(1,5)]]
    end=[[Fraction(2,3)],[Fraction(3,4)]]
    marked=[row[:] for row in bare]; marked[0][0]*=Fraction(2,5); marked[0][1]*=Fraction(7,3)
    exact=sum(marked[0][j]*end[j][0] for j in range(2))
    actual=mm(tolog(np.asarray(marked,dtype=float)),tolog(np.asarray(end,dtype=float)))[0,0]
    assert abs(actual-math.log(float(exact)))<1e-14
    # alpha1 tile grouping preserves the full finite marked product.
    logarithm=tolog(ma); last=tolog(ta)
    full,tail=double_pair(logarithm,last,3)
    for split in (0,1,2):
        a,b=double_pair(logarithm,last,split)
        a,b=double_pair(a,b,3-split)
        assert np.allclose(a,full,atol=1e-14) and np.allclose(b,tail,atol=1e-14)
    floors=0
    for block_sizes in itertools.product(range(5),repeat=4):
        for occupancy in ((0,0,0,0),(1,1,1,1),(0,2,5,16),(16,16,16,16)):
            actual=sum(max(0,13-l) for h,l in zip(block_sizes,occupancy) if h)
            assert actual>=13*((sum(block_sizes)+3)//4)-sum(occupancy)
            floors+=1
    assert (LIMIT+256*256)//13==61503
    return {"exact_route_mask_entries":12,"Jensen_inequalities":12,
            "honest_terminal_mark_tests":1,"alpha1_tile_invariance":3,
            "deterministic_global_floor_checks":floors}


def optimize(z,alpha,branch,matrix,terminal,degrees,coeff,minimize):
    logd=np.log(np.asarray([math.comb(16,j) for j in range(17)],dtype=float))
    profiles=((128,0),(127,2)); bases=[math.log(128/3968),math.log(127/3967),math.log(2/3967)]
    normalizers=[math.log(math.comb(4096,f)*math.comb(4096-f,s)) for f,s in profiles]
    evaluations=0
    def evaluate(values,detail=False):
        nonlocal evaluations
        guard(); evaluations+=1
        logmark=float(values[0]); logfloor=float(values[4]); regions=[]; tails=[]; witnesses=[]
        for i,(f,s) in enumerate(profiles):
            logu=float(values[1+i]); logw=-math.inf if i==0 else float(values[3])
            u=math.exp(logu); w=0. if i==0 else math.exp(logw); logU=math.log1p(u+w)
            allowed=degrees[:,1]==0 if i==0 else np.ones(len(degrees),dtype=bool)
            logs=np.log(coeff[allowed])+degrees[allowed,0]*logu-32*logU
            if i:
                logs+=degrees[allowed,1]*logw
            weights=np.exp(logs)
            raw=np.tensordot(weights,matrix[allowed],axes=(0,0))
            end=np.tensordot(weights,terminal[allowed],axes=(0,0))
            lm=tolog(raw); lt=tolog(end)
            if branch=="births_ge_513":
                lm[0,1:]+=logmark
            else:
                lm[0,0]+=logmark
            charge=np.asarray([(h+3)//4 for h in range(17)],dtype=float)[:,None]*logfloor
            lm+=charge; lt+=charge
            regional,final=region(lm,lt,alpha,logd)
            regions.append(regional); tails.append(final)
            h=4096*logU-f*logu-(s*logw if s else 0)-normalizers[i]
            witnesses.append({"full":f,"singles":s,"u":u,"w":w,"H_bits":h/LN2,
                "regular00_bits":float(regional[0,0]/LN2),
                "regular16to16_bits":float(regional[16,16]/LN2),
                "terminal0_bits":float(final[0,0]/LN2)})
        path,orientation=chain(regions,tails)
        source=alpha*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
        output=-alpha*LIMIT*math.log(z)
        threshold=513 if branch=="births_ge_513" else 24575
        correction=-alpha*threshold*logmark
        floor_correction=-alpha*61503*logfloor
        routing=128*sum(x["H_bits"] for x in witnesses)*LN2
        source_count=math.log(4096)
        score=(source+output+correction+floor_correction+routing+source_count+path)/LN2
        if not detail:
            return score
        return {"z":z,"alpha":alpha,"branch":branch,"log2_bound":score,"raw_branch_mark":math.exp(logmark),
            "raw_global_floor_mark":math.exp(logfloor),
            "components_bits":{"source":source/LN2,"output_Chernoff":output/LN2,
                "branch_correction":correction/LN2,"global_floor_correction":floor_correction/LN2,
                "routing_normalizer":routing/LN2,
                "finite_two_orientation_sum":path/LN2,"sets_per_orientation":source_count/LN2},
            "orientation_chain_bits":[x/LN2 for x in orientation],"regional_witnesses":witnesses}
    initial=[1. if branch=="births_ge_513" else -1.,*bases,-.1]
    bounds=[(0.,20.) if branch=="births_ge_513" else (-20.,0.),
            (bases[0]-6,bases[0]+6),(bases[1]-6,bases[1]+6),(bases[2]-6,bases[2]+6),(-20.,0.)]
    result=minimize(evaluate,np.asarray(initial),method="L-BFGS-B",bounds=bounds,
        options={"maxiter":12,"maxfun":110,"ftol":1e-10,"gtol":1e-5,"maxls":10})
    final=evaluate(result.x,True)
    final.update({"evaluations":evaluations,"iterations":int(result.nit),
                  "optimizer_success":bool(result.success),"optimizer_message":str(result.message)})
    return final


def main():
    global DEADLINE,np,logsumexp
    if len(sys.argv)!=1:
        raise ValueError("fixed tile16 pilot accepts no extra scans")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); DEADLINE=started+48
    import numpy as numpy
    from scipy.special import logsumexp as lse
    from scipy.optimize import minimize
    np=numpy; logsumexp=lse
    singleton,filtered=imports(); checked=tests()
    shapes,index,degrees,coeff,records,checks=geometry()
    results=[]; refined=0; maximum_array_bytes=0
    for z in (.8,.85):
        matrix,terminal,refined=raw_tables(shapes,index,records,z,singleton,filtered)
        maximum_array_bytes=max(maximum_array_bytes,matrix.nbytes+terminal.nbytes)
        for alpha in (.05,.1):
            for branch in ("quiet_le_24575","births_ge_513"):
                results.append(optimize(z,alpha,branch,matrix,terminal,degrees,coeff,minimize))
    print(json.dumps({"scope":"FLOAT actual consecutive256-group supports; all state paths; selected complementary branches",
        "whole_code_certificate":False,"directed":False,"tile_packets":16,"source_sets":8192,
        "matching_histograms":[[128,0],[127,2]],"charged_regular_steps":32767,
        "results":results,"best_per_branch":[min((r for r in results if r["branch"]==branch),key=lambda x:x["log2_bound"])
            for branch in ("quiet_le_24575","births_ge_513")],
        "tests":{**checked,"exact_packet_route_normalizers":checks},"sorted_shapes":len(shapes),
        "route_records":len(records),"pair_only_singleton_refinements":refined,
        "maximum_local_array_bytes":maximum_array_bytes,"estimated_peak_bytes":128*(1<<20),
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"grouped_pin":PIN,
        "seconds":time.monotonic()-started,"generated_bytes":0},indent=2))


if __name__=="__main__":
    main()
