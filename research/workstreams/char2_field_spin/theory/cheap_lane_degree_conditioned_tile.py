"""A/B test: retain tile route degrees before taking the fractional power.

This fixed FLOAT diagnostic uses consecutive256-group supports with actual
alternating histograms(128,0),(127,2), including both orientations. Raw
caps, source count, terminal, tile16, output z=.85, alpha=.05, and quiet
mark .0044813232524843144 match the pinned marked-state pilot. The branch
is quiet regular steps<=24575. The global-floor mark is fixed to1.

Let C[f,s] be the raw16-packet matrix coefficient with f full pairs and
s singleton pairs; it already sums all interior masks and includes final
mask multiplicity D_j=C(16,j). Route normalization is the multinomial
V[f,s]=512!/(f!s!(512-f-s)!). Coefficientwise Jensen gives
 J[f,s,h,j]=V[f,s]*D_j*(C[f,s,h,j]/(V[f,s]*D_j))**alpha.
The terminal coefficient of P**15*T uses the same V but omits D. The
output-only T carries no quiet/birth/feedback factor. Every raw quiet00
factor is included before polynomial multiplication and fractional power.

Raw positive log-matrix convolutions retain f<=128,s<=2; discarded degrees
cannot contribute to either target regional coefficient. No FFT or global
scale is used. Each output matrix entry has its own stable log-sum-exp.
Only after these tile arrays are frozen are route markers optimized.
At markers u,w, evaluate J and divide by(1+u+w)**512. Eight tiles form a
region. The regional Cauchy factor H and actual alternating finite chain
are exactly those of the pinned pilot. This normalization is unchanged.

Compare evaluated-route Jensen and degree-conditioned Jensen at the SAME
old witnesses, then at the SAME newly optimized witnesses. No marks other
than the three route markers are optimized. Improvement combines degree
conditioning with safe impossible-degree truncation inside each tile.
A positive result diagnoses a
relaxation; it does not establish construction failure or exhaust parameters.
FLOAT only, one48-second invocation,20GiB free, estimated peak<128MiB.
No frozen edits, output files, caches, benchmarks, scans, or automatic retry.
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
from types import SimpleNamespace

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
PIN="7ca2f572886d3b19c6a37ed2513cde9c29895e9c8424219f8f547317db0d09b8"
Z=.85
ALPHA=.05
QUIET=.0044813232524843144
OLD=(.032324150838384455,.03193799480012414,.00036989925381939663)
LN2=math.log(2)
DEADLINE=0.
np=None
lse=None


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("48-second deadline exhausted; stop without retry")


def load():
    path=Path(__file__).with_name("cheap_lane_marked_state_pilot.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen marked-state pilot hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_degree_marked",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.np=np; module.logsumexp=lse; module.DEADLINE=DEADLINE
    return module


def convolution(a,b,limit_f=128,limit_s=2):
    """Positive matrix polynomial convolution, per-entry log scales."""
    last_s=min(limit_s,a.shape[0]+b.shape[0]-2)
    last_f=min(limit_f,a.shape[1]+b.shape[1]-2)
    out=np.full((last_s+1,last_f+1,a.shape[2],b.shape[3]),-math.inf)
    for s in range(last_s+1):
        for f in range(last_f+1):
            guard()
            lo=max(0,f-b.shape[1]+1); hi=min(f,a.shape[1]-1)
            ids=np.arange(lo,hi+1); other=f-ids
            for sa in range(max(0,s-b.shape[0]+1),min(s,a.shape[0]-1)+1):
                # Dimensions: split full count, incoming, middle, outgoing.
                terms=a[sa,ids,:,:,None]+b[s-sa,other,None,:,:]
                value=lse(terms,axis=(0,2))
                out[s,f]=np.logaddexp(out[s,f],value)
    assert not np.isnan(out).any()
    return out


def raw_polynomials(m,t,degrees,coeff,pilot):
    p=np.zeros((3,33,17,17)); end=np.zeros((3,33,17,1))
    full=degrees[:,0].astype(np.int16); singles=degrees[:,1].astype(np.int16)
    np.add.at(p,(singles,full),coeff[:,None,None]*m)
    np.add.at(end,(singles,full),coeff[:,None,None]*t)
    p=pilot.tolog(p); end=pilot.tolog(end)
    p[:,:,0,0]+=math.log(QUIET)
    return p,end


def normalize(poly,slots,alpha,multiplicities):
    out=poly.copy(); logd=np.log(np.asarray(multiplicities,dtype=float))
    for s in range(len(poly)):
        for f in range(poly.shape[1]):
            if f+s>slots:
                out[s,f]=-math.inf
            else:
                logv=math.log(math.comb(slots,f)*math.comb(slots-f,s))
                out[s,f]=alpha*poly[s,f]+(1-alpha)*(logv+logd[None,:])
    return out


def exact_mul(a,b,limit_f,limit_s):
    out={}; rows=len(next(iter(a.values()))); columns=len(next(iter(b.values()))[0])
    for (fa,sa),ma in a.items():
        for (fb,sb),mb in b.items():
            if fa+fb>limit_f or sa+sb>limit_s:
                continue
            key=(fa+fb,sa+sb)
            result=out.setdefault(key,[[Fraction(0) for _ in range(columns)] for _ in range(rows)])
            for i,k,j in itertools.product(range(rows),range(len(mb)),range(columns)):
                result[i][j]+=ma[i][k]*mb[k][j]
    return out


def to_logs(poly,pilot):
    first=next(iter(poly.values())); rows=len(first); columns=len(first[0])
    values=np.zeros((max(s for _,s in poly)+1,max(f for f,_ in poly)+1,rows,columns))
    for (f,s),matrix in poly.items():
        values[s,f]=np.asarray(matrix,dtype=float)
    return pilot.tolog(values)


def tests(pilot):
    # One toy slot with empty/full/single choices, two state categories.
    d=(1,2); p={}; t={}
    for f,s in ((0,0),(1,0),(0,1)):
        raw=[[Fraction(1+i+j+f+s,9+2*i+3*j+f+s)*d[j] for j in range(2)] for i in range(2)]
        raw[0][0]*=Fraction(3,7)
        p[(f,s)]=raw
        t[(f,s)]=[[Fraction(1+i+f+s,11+i+f+s)] for i in range(2)]
    a,b=to_logs(p,pilot),to_logs(t,pilot); exact_a,exact_b=p,t; comparisons=0
    for packets in (2,4):
        exact_b,exact_a=exact_mul(exact_a,exact_b,4,2),exact_mul(exact_a,exact_a,4,2)
        b,a=convolution(a,b,4,2),convolution(a,a,4,2)
        for exact,actual in ((exact_a,a),(exact_b,b)):
            for (f,s),matrix in exact.items():
                expected=np.log(np.asarray(matrix,dtype=float))
                assert np.allclose(expected,actual[s,f],atol=4e-13,rtol=4e-13)
                comparisons+=expected.size
        assert np.array_equal(normalize(a,packets,1.,d),a)
        assert np.array_equal(normalize(b,packets,1.,[1]),b)
    # Positive Jensen, with route masses retained outside alpha.
    u,w=Fraction(1,3),Fraction(1,5); route_mass=(1+u+w)**4; alpha=.3
    powered=normalize(a,4,alpha,d); powered_t=normalize(b,4,alpha,[1])
    routes=np.asarray([[float(u**f*w**s) for f in range(a.shape[1])] for s in range(a.shape[0])])
    raw=(np.exp(a)*routes[:,:,None,None]).sum(axis=(0,1))
    conditioned=(np.exp(powered)*routes[:,:,None,None]).sum(axis=(0,1))
    upper=raw**alpha*(float(route_mass)*np.asarray(d)[None,:])**(1-alpha)
    assert np.all(conditioned<=upper*(1+2e-13))
    rawt=(np.exp(b)*routes[:,:,None,None]).sum(axis=(0,1))
    conditioned_t=(np.exp(powered_t)*routes[:,:,None,None]).sum(axis=(0,1))
    assert np.all(conditioned_t<=rawt**alpha*float(route_mass)**(1-alpha)*(1+2e-13))
    return {"exact_bivariate_matrix_entries":comparisons,"alpha1_invariance":4,
            "degree_vs_evaluated_Jensen_inequalities":6,"terminal_matrix_order":"passed"}


def tile_value(poly,logu,logw,logU):
    singles=range(len(poly)) if math.isfinite(logw) else range(1)
    terms=[]
    full=np.arange(poly.shape[1],dtype=float)
    for s in singles:
        terms.append(poly[s]+(full*logu+(s*logw if s else 0))[:,None,None])
    return lse(np.concatenate(terms,axis=0),axis=0)-512*logU


def evaluator(m,t,degrees,coeff,a,b,pilot):
    logd=np.log(np.asarray([math.comb(16,j) for j in range(17)],dtype=float))
    profiles=((128,0),(127,2))
    norms=[math.log(math.comb(4096,f)*math.comb(4096-f,s)) for f,s in profiles]
    source=ALPHA*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
    output=-ALPHA*pilot.LIMIT*math.log(Z)
    correction=-ALPHA*24575*math.log(QUIET)
    def evaluate(values,kind="conditioned",detail=False):
        guard(); regions=[]; tails=[]; witnesses=[]
        for i,(f,s) in enumerate(profiles):
            logu=float(values[i]); logw=-math.inf if i==0 else float(values[2])
            u=math.exp(logu); w=0. if i==0 else math.exp(logw); logU=math.log1p(u+w)
            if kind=="conditioned":
                ta=tile_value(a,logu,logw,logU); tb=tile_value(b,logu,logw,logU)
                reg,end=pilot.double_pair(ta,tb,3)
            else:
                allowed=degrees[:,1]==0 if i==0 else np.ones(len(degrees),dtype=bool)
                logs=np.log(coeff[allowed])+degrees[allowed,0]*logu-32*logU
                if i:
                    logs+=degrees[allowed,1]*logw
                weights=np.exp(logs)
                raw=pilot.tolog(np.tensordot(weights,m[allowed],axes=(0,0)))
                terminal=pilot.tolog(np.tensordot(weights,t[allowed],axes=(0,0)))
                raw[0,0]+=math.log(QUIET)
                reg,end=pilot.region(raw,terminal,ALPHA,logd)
            regions.append(reg); tails.append(end)
            h=4096*logU-f*logu-(s*logw if s else 0)-norms[i]
            witnesses.append({"full":f,"single":s,"u":u,"w":w,"H_bits":h/LN2,
                "regular00_bits":float(reg[0,0]/LN2),"regular16to16_bits":float(reg[16,16]/LN2),
                "terminal0_bits":float(end[0,0]/LN2)})
        paths,orientations=pilot.chain(regions,tails)
        routing=128*sum(row["H_bits"] for row in witnesses)*LN2
        score=(source+output+correction+routing+math.log(4096)+paths)/LN2
        if not detail:
            return score
        return {"kind":kind,"log2_bound":score,"regional_witnesses":witnesses,
            "orientation_chain_bits":[v/LN2 for v in orientations],
            "components_bits":{"source":source/LN2,"output_Chernoff":output/LN2,
                "quiet_correction":correction/LN2,"routing_normalizer":routing/LN2,
                "sets_per_orientation":12.,"finite_two_orientation_sum":paths/LN2}}
    return evaluate


def main():
    global DEADLINE,np,lse
    if len(sys.argv)!=1:
        raise ValueError("fixed A/B test accepts no scans")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); DEADLINE=started+48
    import numpy as numpy
    from scipy.special import logsumexp
    from scipy.optimize import minimize
    np=numpy; lse=logsumexp
    pilot=load(); checked=tests(pilot)
    singleton,filtered=pilot.imports()
    shapes,index,degrees,coeff,records,count_checks=pilot.geometry()
    m,t,refined=pilot.raw_tables(shapes,index,records,Z,singleton,filtered)
    a,b=raw_polynomials(m,t,degrees,coeff,pilot); phases=[]
    for stage in range(4):
        start=time.monotonic()
        b,a=convolution(a,b),convolution(a,a)
        phases.append({"packets":2**(stage+1),"seconds":time.monotonic()-start,
                       "finite_matrix_coefficients":int(np.isfinite(a).sum())})
        print(json.dumps({"progress":"raw_tile_complete",**phases[-1],
                          "total_seconds":time.monotonic()-started}),flush=True)
    raw_array_bytes=a.nbytes+b.nbytes
    a=normalize(a,512,ALPHA,[math.comb(16,j) for j in range(17)])
    b=normalize(b,512,ALPHA,[1])
    evaluate=evaluator(m,t,degrees,coeff,a,b,pilot)
    old_logs=np.log(np.asarray(OLD)); old_ab={kind:evaluate(old_logs,kind,True) for kind in ("evaluated","conditioned")}
    assert abs(old_ab["evaluated"]["log2_bound"]-37827.01773005103)<1e-6
    assert old_ab["conditioned"]["log2_bound"]<=old_ab["evaluated"]["log2_bound"]+1e-8
    print(json.dumps({"progress":"fixed_witness_AB","evaluated_log2":old_ab["evaluated"]["log2_bound"],
        "conditioned_log2":old_ab["conditioned"]["log2_bound"],"seconds":time.monotonic()-started}),flush=True)
    calls=0
    def objective(values):
        nonlocal calls
        calls+=1
        return evaluate(values)
    if DEADLINE-time.monotonic()>=8:
        optimized=minimize(objective,old_logs,method="L-BFGS-B",bounds=[(x-5,x+5) for x in old_logs],
                           options={"maxiter":20,"maxfun":120,"ftol":1e-11,"gtol":1e-6,"maxls":12})
    else:
        optimized=SimpleNamespace(x=old_logs,success=False,nit=0,message="optimization skipped to preserve deadline")
    new_ab={kind:evaluate(optimized.x,kind,True) for kind in ("evaluated","conditioned")}
    assert new_ab["conditioned"]["log2_bound"]<=new_ab["evaluated"]["log2_bound"]+1e-8
    print(json.dumps({"scope":"FLOAT actual consecutive256-group support, quiet<=24575; degree-conditioned tile16 A/B",
        "whole_code_certificate":False,"directed":False,"parameters":{"z":Z,"alpha":ALPHA,"quiet_mark":QUIET,"floor_mark":1},
        "same_prior_witness":old_ab,"same_new_witness":new_ab,"raw_convolution_stages":phases,
        "optimizer":{"success":bool(optimized.success),"message":str(optimized.message),"calls":calls,"iterations":int(optimized.nit)},
        "tests":{**checked,"exact_packet_route_normalizers":count_checks},"pair_only_singleton_refinements":refined,
        "raw_tile_bytes":raw_array_bytes,"estimated_peak_bytes":128*(1<<20),"charged_regular_steps":32767,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"marked_pilot_pin":PIN,
        "seconds":time.monotonic()-started,"generated_bytes":0},indent=2))


if __name__=="__main__":
    main()
