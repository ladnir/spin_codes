"""Fixed FLOAT A/B: retain actual active RS-block count in the state class.

This tests only consecutive256-group supports, with actual alternating
profiles(128full,0single)/(127full,2single). It is not a certificate or a
parameter search. z=.85,alpha=.05,tile16,quiet mark and route witnesses are
copied from the hash-pinned17-state pilot. The branch is Q<=24575; its
birth restriction is dropped. The global-floor mark is1.

States(h,B) record actual F-state support size h and number B of active
four-coordinate RS blocks. They do NOT assert auxiliary E-state support.
D(h,B)=C(4,B)[x^h]((1+x)^4-1)^B. There are35 classes, including(0,0).
For each fixed occupancy, maximize only over t in{0,..,4}^4 with sum h
and exactly B positive entries. No C(4,t_i) appears in this MAX. Returns
require t_i=0 at empty input blocks. Other transitions use the necessary
support condition j+4*occupied+sum_{ell_i=0}t_i>=17.

Reuse the frozen filtered scalar block moments, gamma choices and joint
live filter. For a fixed t, retain the actual floor sum_active(13-ell_i)+,
feedback image lower rank k=sum_i(min(4,ell_i)-t_i)+, and outside-output
rank a=sum_i(t_i-(ell_i-12)+)+. Minimize the compatible scalar caps before
maximizing over t. Then cap every entry with the OLD uniform-h kernel.
The old shared-point singleton extinction bound is consequently retained
without replacing joint point averages by independent-block averages.

The per-fixed-destination cap R[(h,B),j] ignores destination block count.
With L[j,(j,C)]=D(j,C), the full raw matrix is R L. Products factor through
17 dimensions: (R L)^p=R(L R)^(p-1)L. Endpoint-normalized Jensen preserves
this form: after alpha, the endpoint D(j,C) remains outside the power.
Terminal R(LR)^14 L T contains15 updates and one output-only packet.
There is no final destination multiplicity, quiet mark or feedback.

Also report the0/full16-only channel, restricting raw states BEFORE tile
pooling. Those classes must agree with the old kernel. All finite chains
include32767 updates, both source orientations, and their positive sum.
No optimization, scans, caches or result files.48seconds,20GiB free,
estimated peak<128MiB. Implementation may be inspected before authorization
to run; importing this file performs no computation.
"""

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
PIN="7ca2f572886d3b19c6a37ed2513cde9c29895e9c8424219f8f547317db0d09b8"
Z=.85
ALPHA=.05
QUIET=.0044813232524843144
ROUTE=(.032324150838384455,.03193799480012414,.00036989925381939663)
Q=256
LN2=math.log(2)
DEADLINE=0.
np=None
lse=None


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("48-second bound exhausted; stop without retry")


def load():
    path=Path(__file__).with_name("cheap_lane_marked_state_pilot.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen marked-state pilot hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_active_blocks_pilot",path)
    pilot=importlib.util.module_from_spec(spec); spec.loader.exec_module(pilot)
    pilot.np=np; pilot.logsumexp=lse; pilot.DEADLINE=DEADLINE
    return pilot


def classes():
    states=[(0,0)]+[(h,b) for h in range(1,17) for b in range(1,5) if b<=h<=4*b]
    index={state:i for i,state in enumerate(states)}
    counts=[0]*len(states)
    placements=[]
    for t in itertools.product(range(5),repeat=4):
        h=sum(t); b=sum(x>0 for x in t)
        counts[index[h,b]]+=math.prod(math.comb(4,x) for x in t)
        placements.append((t,index[h,b],h,b))
    assert len(states)==35 and sum(counts)==65536
    for h in range(17):
        assert sum(counts[i] for i,state in enumerate(states) if state[0]==h)==math.comb(16,h)
    assert [counts[index[8,b]] for b in (2,3,4)]==[6,1968,10896]
    return states,index,counts,placements


def block_tables(residual,floor):
    """Same retained-Z caps as the frozen scalar helper, before alpha."""
    g=residual+(1-residual)/Q
    out=np.zeros((17,5)); keep=np.ones((17,5))
    for ell in range(17):
        n=16-ell
        weights=[math.comb(n,v)*residual**(n-v)*(1-residual)**v
                 for v in range(min(3,n)+1)]
        retained=sum(weights)
        out[ell,0]=g**ell
        for h in range(1,5):
            rank=sum(w*Q**(-max(0,h-4+v)) for v,w in enumerate(weights))
            split=Q**(-h)*retained+(1-Q**(-h))*residual**max(0,n-3)
            value=g**ell*min(rank,split)
            if h==1:
                value=min(value,g**ell*residual**max(0,n-2))
            factor=floor**max(0,13-ell)
            out[ell,h]=factor*value
            keep[ell,h]=factor*retained
    return out,keep


def refined_tables(occ,states,placements):
    """A class-uniform raw cap: max_t min(valid compatible scalar caps)."""
    n=len(occ); ell=occ.sum(axis=1); ci=np.minimum(occ,4)
    c=ci.sum(axis=1); occupied=(occ>0).sum(axis=1)
    dest=np.arange(17); missing=16-dest
    joint=np.asarray([1.]+[min(1.,Q/(Q-1)*255.**(-(16-j))) for j in range(1,17)])
    result=terminal=None
    for gamma in (0.,.5,1.):
        floor=Z**gamma; residual=Z**(1-gamma)
        g=residual+(1-residual)/Q
        bout,bkeep=block_tables(residual,floor)
        current=np.zeros((n,35,17)); end=np.zeros((n,35,1))
        # These rows are taken unchanged from the old exact quiet/activation cap.
        current[:,0,:]=math.inf; end[:,0,0]=math.inf
        for raw_t,cell,h,b in placements:
            guard()
            if h==0:
                continue
            t=np.asarray(raw_t,dtype=np.int16)
            output=np.prod(bout[occ,t[None,:]],axis=1)
            retained=np.prod(bkeep[occ,t[None,:]],axis=1)
            local_floor=np.where(t[None,:]>0,np.maximum(0,13-occ),0).sum(axis=1)
            k=np.maximum(0,ci-t[None,:]).sum(axis=1)
            outside_rank=np.maximum(0,t[None,:]-np.maximum(0,occ-12)).sum(axis=1)
            rf=np.minimum(missing,h)
            d0=np.maximum(0,c[:,None]-dest[None,:]-rf[None,:])
            e0=np.maximum(0,outside_rank[:,None]-rf[None,:])
            padding=(np.power(floor,local_floor)[:,None]*np.power(float(Q),-rf)[None,:]
                *np.power(g,ell)[:,None]*np.power(Q*g,-d0)*np.power(g,e0))
            # Restrict a fixed k-dimensional C-image before averaging D1.
            qk=np.power(float(Q),-k)
            fm=np.ones((n,17)); fm[:,0]=qk
            for j in range(1,17):
                if h<missing[j]:
                    fm[:,j]=qk+(1-qk)*255.**(-(missing[j]-h))
            feedback=np.power(float(Q),-rf)[None,:]*fm
            value=np.minimum(padding,np.minimum(output[:,None],feedback*retained[:,None]))
            value[:,1:]=np.minimum(value[:,1:],output[:,None]*joint[None,1:])
            union_size=4*occupied+np.where(occ==0,t[None,:],0).sum(axis=1)
            feasible=union_size[:,None]+dest[None,:]>=17
            feasible[:,0]=np.all((t[None,:]==0)|(occ>0),axis=1)
            value*=feasible
            current[:,cell,:]=np.maximum(current[:,cell,:],value)
            end[:,cell,0]=np.maximum(end[:,cell,0],output)
        result=current if result is None else np.minimum(result,current)
        terminal=end if terminal is None else np.minimum(terminal,end)
    assert np.isfinite(result[:,1:]).all() and (result>=0).all()
    assert np.isfinite(terminal[:,1:]).all() and (terminal>=0).all()
    return result,terminal


def exact_mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def tests(pilot,states,index,counts,placements):
    # Positive block-count enumeration agrees with the displayed polynomial.
    for b in range(1,5):
        polynomial=[1]
        for _ in range(b):
            out=[0]*(len(polynomial)+4)
            for i,x in enumerate(polynomial):
                for j in range(1,5): out[i+j]+=x*math.comb(4,j)
            polynomial=out
        for h in range(b,4*b+1):
            assert counts[index[h,b]]==math.comb(4,b)*polynomial[h]
    checks=0
    for raw_t,_,h,b in placements:
        for occ in ((0,0,0,0),(1,0,0,0),(1,1,0,0),(0,2,5,16),(16,16,16,16)):
            c=sum(min(4,x) for x in occ)
            k=sum(max(0,min(4,x)-t) for x,t in zip(occ,raw_t))
            a=sum(max(0,t-max(0,x-12)) for x,t in zip(occ,raw_t))
            old_a=max(0,h-sum(max(0,x-12) for x in occ))
            floor=sum(max(0,13-x) for x,t in zip(occ,raw_t) if t)
            assert max(0,c-h)<=k<=16-h and old_a<=a<=h
            assert floor>=13*b-sum(occ)
            if all(t==0 or x>0 for x,t in zip(occ,raw_t)):
                assert b<=sum(x>0 for x in occ) and h<=4*sum(x>0 for x in occ)
            checks+=1
    # Exact toy factorization has two classes with the same cardinality.
    hs=(0,1,2,2); ds=(1,2,1,3)
    r=[[Fraction(1+i+j,9+2*i+3*j) for j in range(3)] for i in range(4)]
    l=[[Fraction(ds[j] if hs[j]==i else 0) for j in range(4)] for i in range(3)]
    t=[[Fraction(2+i,11+i)] for i in range(4)]
    r[0][0]*=Fraction(3,7)  # Only the regular update gets this quiet mark.
    m=exact_mm(r,l); lr=exact_mm(l,r)
    raw=exact_mm(exact_mm(m,m),m)
    fixed=exact_mm(exact_mm(r,lr),lr)
    assert raw==exact_mm(fixed,l)
    end=exact_mm(exact_mm(m,m),t)
    assert end==exact_mm(exact_mm(r,lr),exact_mm(l,t))
    alpha=.3
    full=np.asarray(raw,dtype=float)**alpha*np.asarray(ds)[None,:]**(1-alpha)
    fact=np.asarray(fixed,dtype=float)**alpha@np.asarray(l,dtype=float)
    assert np.allclose(full,fact,rtol=3e-15,atol=0)
    log_r=pilot.tolog(np.asarray(r,dtype=float)); log_l=pilot.tolog(np.asarray(l,dtype=float))
    log_t=pilot.tolog(np.asarray(t,dtype=float)); log_lr=pilot.mm(log_l,log_r)
    terminal=pilot.mm(pilot.mm(log_r,log_lr),pilot.mm(log_l,log_t))
    assert np.allclose(terminal,np.log(np.asarray(end,dtype=float)),atol=2e-14,rtol=0)
    assert 255*128+127==32767
    return {"exact_class_counts":35,"per_h_multiplicity_identity":17,
        "placement_rank_floor_checks":checks,"exact_factorized_raw_entries":16,
        "endpoint_Jensen_factorization":16,"honest_terminal_entries":4,
        "regular_update_count":32767}


def factored_region(raw,end,states,counts,pilot):
    """35 classes,17-dimensional raw products, honest terminal tile."""
    h=np.asarray([x[0] for x in states],dtype=int)
    log_l=np.full((17,35),-math.inf)
    log_l[h,np.arange(35)]=np.log(np.asarray(counts,dtype=float))
    log_r=pilot.tolog(raw); log_t=pilot.tolog(end)
    log_r[0,0]+=math.log(QUIET)
    lr=pilot.mm(log_l,log_r)
    tile_r=ALPHA*pilot.mm(log_r,pilot.power(lr,15))
    tile_t=ALPHA*pilot.mm(pilot.mm(log_r,pilot.power(lr,14)),pilot.mm(log_l,log_t))
    pooled=pilot.mm(log_l,tile_r)
    regional_r=pilot.mm(tile_r,pilot.power(pooled,7))
    regional_t=pilot.mm(pilot.mm(tile_r,pilot.power(pooled,6)),pilot.mm(log_l,tile_t))
    return pilot.mm(regional_r,log_l),regional_t


def evaluate(pilot,old_m,old_t,refined,refined_t,shape_index,degrees,coeff,states,counts):
    profiles=((128,0),(127,2)); d=np.asarray([math.comb(16,j) for j in range(17)],dtype=float)
    all_results={name:([],[]) for name in ("old17","refined35","old_0full16","refined_0full16")}
    witnesses=[]; full=states.index((16,4))
    for matching,(f,s) in enumerate(profiles):
        guard()
        u=ROUTE[matching]; w=0. if matching==0 else ROUTE[2]; logU=math.log1p(u+w)
        allowed=degrees[:,1]==0 if matching==0 else np.ones(len(degrees),dtype=bool)
        logs=np.log(coeff[allowed])+degrees[allowed,0]*math.log(u)-32*logU
        if matching: logs+=degrees[allowed,1]*math.log(w)
        weights=np.exp(logs)
        raw_old=np.tensordot(weights,old_m[allowed],axes=(0,0))
        end_old=np.tensordot(weights,old_t[allowed],axes=(0,0))
        raw_new=np.empty((35,17)); end_new=np.empty((35,1))
        rows=shape_index[allowed]
        for i,(h,_) in enumerate(states):
            # Each record retains the frozen shared-point singleton bound.
            cap=np.minimum(refined[rows,i,:],old_m[allowed,h,:]/d[None,:])
            terminal=np.minimum(refined_t[rows,i,0],old_t[allowed,h,0])
            raw_new[i]=weights@cap; end_new[i,0]=weights@terminal
            assert np.all(raw_new[i]<=raw_old[h]/d*(1+2e-14))
        # These extreme classes have no class split and must remain unchanged.
        for old_h,new_h in ((0,0),(16,full)):
            assert np.allclose(raw_new[new_h],raw_old[old_h]/d,rtol=3e-13,atol=1e-300)
            assert math.isclose(float(end_new[new_h,0]),float(end_old[old_h,0]),rel_tol=3e-13)
        lm=pilot.tolog(raw_old); lm[0,0]+=math.log(QUIET)
        reg,tail=pilot.region(lm,pilot.tolog(end_old),ALPHA,np.log(d))
        all_results["old17"][0].append(reg); all_results["old17"][1].append(tail)
        reg,tail=factored_region(raw_new,end_new,states,counts,pilot)
        all_results["refined35"][0].append(reg); all_results["refined35"][1].append(tail)
        for name,raw,end in (("old_0full16",raw_old[np.ix_([0,16],[0,16])],end_old[[0,16]]),
                ("refined_0full16",raw_new[np.ix_([0,full],[0,16])],end_new[[0,full]])):
            lm=pilot.tolog(raw); lm[0,0]+=math.log(QUIET)
            reg,tail=pilot.region(lm,pilot.tolog(end),ALPHA,np.zeros(2))
            all_results[name][0].append(reg); all_results[name][1].append(tail)
        normalizer=math.log(math.comb(4096,f)*math.comb(4096-f,s))
        h=4096*logU-f*math.log(u)-(s*math.log(w) if s else 0)-normalizer
        witnesses.append({"full":f,"single":s,"u":u,"w":w,"H_bits":h/LN2})
    source=ALPHA*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
    output=-ALPHA*pilot.LIMIT*math.log(Z)
    correction=-ALPHA*24575*math.log(QUIET)
    routing=128*sum(x["H_bits"] for x in witnesses)*LN2
    prefix=source+output+correction+routing+math.log(4096)
    results={}
    for name,(regions,tails) in all_results.items():
        path,orientations=pilot.chain(regions,tails)
        results[name]={"log2_bound":(prefix+path)/LN2,
            "finite_two_orientation_sum_bits":path/LN2,
            "orientation_chain_bits":[x/LN2 for x in orientations],
            "regional00_bits":[float(x[0,0]/LN2) for x in regions]}
    assert abs(results["old17"]["log2_bound"]-37827.01773005103)<1e-6
    assert results["refined35"]["log2_bound"]<=results["old17"]["log2_bound"]+1e-6
    assert abs(results["old_0full16"]["log2_bound"]-results["refined_0full16"]["log2_bound"])<1e-6
    return {"results":results,"refinement_change_bits":results["refined35"]["log2_bound"]-results["old17"]["log2_bound"],
        "fixed_route_witnesses":witnesses,"components_bits":{"source":source/LN2,
        "output_Chernoff":output/LN2,"quiet_mark_correction":correction/LN2,
        "routing_normalizer":routing/LN2,"source_sets_per_orientation":12.}}


def main():
    global DEADLINE,np,lse
    if len(sys.argv)!=1:
        raise ValueError("fixed A/B helper accepts no scans or optimization arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the source volume")
    started=time.monotonic(); DEADLINE=started+48
    import numpy as numpy
    from scipy.special import logsumexp
    np=numpy; lse=logsumexp
    pilot=load(); singleton,filtered=pilot.imports()
    states,index,counts,placements=classes(); checked=tests(pilot,states,index,counts,placements)
    shapes,shape_index,degrees,coeff,records,geometry_checks=pilot.geometry()
    old_m,old_t,singleton_refinements=pilot.raw_tables(shapes,shape_index,records,Z,singleton,filtered)
    estimated_peak=(64*(1<<20)+old_m.nbytes+old_t.nbytes+len(shapes)*35*18*8*3)
    if estimated_peak>=128*(1<<20):
        raise RuntimeError("bounded128MiB estimated peak exceeded before refinement")
    refined,refined_t=refined_tables(shapes,states,placements)
    result=evaluate(pilot,old_m,old_t,refined,refined_t,shape_index,degrees,coeff,states,counts)
    guard()
    print(json.dumps({"scope":"FLOAT fixed35-state A/B on actual consecutive256-group supports",
        "whole_code_certificate":False,"directed":False,"encoder_changed":False,
        "parameters":{"z":Z,"alpha":ALPHA,"tile":16,"quiet_mark":QUIET,"global_floor_mark":1},
        "branch":"quiet regular steps<=24575; birth restriction omitted",
        "states":[{"h":h,"active_blocks":b,"multiplicity":counts[i]} for i,(h,b) in enumerate(states)],
        **result,"tests":{**checked,"exact_packet_route_normalizers":geometry_checks,
            "0_full16_local_identity":"passed","frozen_old_score_replayed":"passed"},
        "sorted_shapes":len(shapes),"route_records":len(records),"singleton_refinements":singleton_refinements,
        "estimated_peak_bytes":estimated_peak,"generated_bytes":0,"seconds":time.monotonic()-started,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"marked_pilot_pin":PIN},indent=2))


if __name__=="__main__":
    main()
