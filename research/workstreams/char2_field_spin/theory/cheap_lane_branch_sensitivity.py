"""Positive dual sensitivity of the exact regional zero/full16 control.

Keep the frozen control's raw caps, z=.85, alpha=.05, quiet mark,
floor1, consecutive256-group profiles, and exact regional coefficients.
The two actual states are zero/full16 throughout. This remains a positive
submajorant, not an upper bound on all state paths or a certificate.

Carry four raw log derivatives: births k on0->16; quiet q on0->0;
active incoming packets a on row16, INCLUDING the output-only terminal;
and the final incoming16 indicator. Regular products add these counts.
Each coefficient stores its log mass and normalized nonnegative derivative
means. Multiplication uses positive log-sum-exp weights, then averages the
sum of the factor means. No finite differences or optimization are used.

At exact regional degree, the normalized coefficient is raised to alpha.
Its log derivatives are multiplied by alpha, then propagated through the
actual finite256-region chain. Divide the final log derivatives by alpha
to report hierarchical tilted means. They satisfy
 E[q]+E[k]+E[a]-Pr(final incoming16)=32767.

For a raw mark exp(x), the corrected score derivative in bits is
 alpha*(mean count-threshold)/log(2). Restore k<=512 with x_k<=0;
retain q<=24575 with its existing negative x_q; and test4a<=61503 with
a floor mark x_floor<=0. The floor slope uses4E[a], not E[a].
At x=0 a negative slope excludes improvement by moving that coordinate
negative while holding the others fixed. It does not exhaust joint changes
to other parameters. Every reported value is FLOAT.

One48-second job,20GiB free,128MiB estimated peak, stdout only. No scans,
new encoder work, frozen edits, generated files, caches, or retries.
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
PIN="e2dbb5fdd206feb0453dee426c17f3eb95283433542e411733c65238f0175586"
LN2=math.log(2)
DEADLINE=0.
np=None
lse=None


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("48-second dual diagnostic exhausted; no retry")


def import_pinned(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen "+label+" hash mismatch")
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def load():
    exact=import_pinned("cheap_lane_exact_region_control.py",PIN,"exact_region_control")
    pool=import_pinned("cheap_lane_pooling_scale_pilot.py",exact.PIN,"dual_pool")
    pool.np=np; pool.lse=lse; pool.DEADLINE=DEADLINE
    degree,pilot=pool.load()
    return pool,degree,pilot


def weighted(logs,means,axes):
    """All normalized derivatives are positive weighted averages."""
    total=lse(logs,axis=axes,keepdims=True)
    weights=np.exp(logs-np.where(np.isfinite(total),total,0.))
    result=(weights[...,None]*means).sum(axis=axes)
    return np.squeeze(total,axis=axes),result


def convolution(a,b,limit_f=128,limit_s=2):
    av,ad=a; bv,bd=b
    last_s=min(limit_s,av.shape[0]+bv.shape[0]-2)
    last_f=min(limit_f,av.shape[1]+bv.shape[1]-2)
    shape=(last_s+1,last_f+1,av.shape[2],bv.shape[3])
    cv=np.full(shape,-math.inf); cd=np.zeros(shape+(4,))
    for s in range(last_s+1):
        for f in range(last_f+1):
            guard(); values=[]; derivatives=[]
            ids=np.arange(max(0,f-bv.shape[1]+1),min(f,av.shape[1]-1)+1)
            other=f-ids
            for sa in range(max(0,s-bv.shape[0]+1),min(s,av.shape[0]-1)+1):
                values.append(av[sa,ids,:,:,None]+bv[s-sa,other,None,:,:])
                derivatives.append(ad[sa,ids,:,:,None,:]+bd[s-sa,other,None,:,:,:])
            cv[s,f],cd[s,f]=weighted(np.concatenate(values,axis=0),np.concatenate(derivatives,axis=0),(0,2))
    assert not np.isnan(cv).any() and np.isfinite(cd).all() and (cd>=0).all()
    return cv,cd


def mm(a,b):
    av,ad=a; bv,bd=b
    return weighted(av[:,:,None]+bv[None,:,:],ad[:,:,None,:]+bd[None,:,:,:],(1,))


def power(a,n):
    rows=a[0].shape[0]
    value=np.full((rows,rows),-math.inf); np.fill_diagonal(value,0.)
    result=(value,np.zeros((rows,rows,4)))
    while n:
        if n&1:
            result=mm(result,a)
        n//=2
        if n:
            a=mm(a,a)
    return result


def chain(regions,tails,pairs=128):
    results=[]
    for i in (0,1):
        value=mm(mm(power(mm(regions[i],regions[1-i]),pairs-1),regions[i]),tails[1-i])
        results.append((float(value[0][0,0]),value[1][0,0]))
    logs=np.asarray([value for value,_ in results]); means=np.asarray([mean for _,mean in results])
    value,mean=weighted(logs,means,(0,))
    return float(value),mean,results


def seed(logs,terminal=False):
    derivatives=np.zeros(logs.shape+(4,))
    if terminal:
        derivatives[:,:,1,:,2]=1.
        derivatives[:,:,1,:,3]=1.
    else:
        derivatives[:,:,0,1,0]=1.
        derivatives[:,:,0,0,1]=1.
        derivatives[:,:,1,:,2]=1.
    return logs,derivatives


def exact_mul(a,b,limit_f,limit_s):
    out={}
    for (fa,sa),(va,da) in a.items():
        for (fb,sb),(vb,db) in b.items():
            if fa+fb>limit_f or sa+sb>limit_s:
                continue
            rows=len(va); columns=len(vb[0])
            vc,dc=out.setdefault((fa+fb,sa+sb),
                ([[Fraction(0) for _ in range(columns)] for _ in range(rows)],
                 [[[Fraction(0) for _ in range(4)] for _ in range(columns)] for _ in range(rows)]))
            for i,k,j in itertools.product(range(rows),range(len(vb)),range(columns)):
                vc[i][j]+=va[i][k]*vb[k][j]
                for axis in range(4):
                    dc[i][j][axis]+=da[i][k][axis]*vb[k][j]+va[i][k]*db[k][j][axis]
    return out


def tests(pilot):
    p={}; t={}
    for f,s in ((0,0),(1,0),(0,1)):
        value=[[Fraction(1+i+j+f+s,9+2*i+3*j+f+s) for j in range(2)] for i in range(2)]
        value[0][0]*=Fraction(3,7)
        terminal=[[Fraction(1+i+f+s,11+i+f+s)] for i in range(2)]
        deriv=[[[value[i][j]*int((axis==0 and i==0 and j==1) or
                    (axis==1 and i==0 and j==0) or (axis==2 and i==1))
                    for axis in range(4)] for j in range(2)] for i in range(2)]
        end_deriv=[[[terminal[i][0]*int(i==1 and axis in (2,3)) for axis in range(4)]] for i in range(2)]
        p[(f,s)]=(value,deriv); t[(f,s)]=(terminal,end_deriv)
    def convert(poly):
        first=next(iter(poly.values()))[0]; rows=len(first); columns=len(first[0])
        values=np.zeros((2,2,rows,columns)); means=np.zeros(values.shape+(4,))
        for (f,s),(v,d) in poly.items():
            values[s,f]=np.asarray(v,dtype=float)
            means[s,f]=np.asarray(d,dtype=float)/values[s,f,:,:,None]
        return pilot.tolog(values),means
    a,b=convert(p),convert(t); ea,eb=p,t; entries=0
    for _ in range(2):
        eb,ea=exact_mul(ea,eb,4,2),exact_mul(ea,ea,4,2)
        b,a=convolution(a,b,4,2),convolution(a,a,4,2)
        for exact,actual in ((ea,a),(eb,b)):
            for (f,s),(v,d) in exact.items():
                va=np.asarray(v,dtype=float); da=np.asarray(d,dtype=float)/va[:,:,None]
                assert np.allclose(actual[0][s,f],np.log(va),atol=4e-13,rtol=4e-13)
                assert np.allclose(actual[1][s,f],da,atol=4e-13,rtol=4e-13)
                entries+=da.size
    # Independent exact sequential product for a four-region finite chain.
    reg0=ea[(1,0)]; reg1=ea[(0,1)]; end0=eb[(1,0)]; end1=eb[(0,1)]
    def matdual(exact):
        v,d=exact; vv=np.asarray(v,dtype=float)
        return np.log(vv),np.asarray(d,dtype=float)/vv[:,:,None]
    value,mean,_=chain([matdual(reg0),matdual(reg1)],[matdual(end0),matdual(end1)],2)
    summed_value=Fraction(0); summed_deriv=[Fraction(0)]*4
    for order in ((reg0,reg1,reg0,end1),(reg1,reg0,reg1,end0)):
        current={(0,0):order[0]}
        for factor in order[1:]:
            current=exact_mul(current,{(0,0):factor},0,0)
        vv,dd=current[(0,0)]; summed_value+=vv[0][0]
        summed_deriv=[a+b for a,b in zip(summed_deriv,dd[0][0])]
    assert abs(value-math.log(float(summed_value)))<1e-12
    assert np.allclose(mean,[float(x/summed_value) for x in summed_deriv],atol=1e-12)
    assert abs(mean[0]+mean[1]+mean[2]-mean[3]-15)<1e-12
    return {"exact_rational_dual_entries":entries,"exact_finite_chain_derivatives":4,
            "exact_finite_chain_count_identity":"passed"}


def main():
    global DEADLINE,np,lse
    if len(sys.argv)!=1:
        raise ValueError("fixed sensitivity diagnostic accepts no scans")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); DEADLINE=started+48
    import numpy as numpy
    from scipy.special import logsumexp
    np=numpy; lse=logsumexp
    pool,degree,pilot=load(); checked=tests(pilot)
    singleton,filtered=pilot.imports()
    shapes,index,degrees,coeff,records,route_checks=pilot.geometry()
    m,t,_=pilot.raw_tables(shapes,index,records,degree.Z,singleton,filtered)
    av,bv=degree.raw_polynomials(m,t,degrees,coeff,pilot)
    av=av[:,:,[0,16],:][:,:,:,[0,16]].copy(); bv=bv[:,:,[0,16],:].copy()
    a,b=seed(av),seed(bv,True)
    del m,t,av,bv
    for _ in range(7):
        b,a=convolution(a,b),convolution(a,a)
    regions=[]; tails=[]
    for f,s in ((128,0),(127,2)):
        norm=math.log(math.comb(4096,f)*math.comb(4096-f,s))
        regions.append((degree.ALPHA*(a[0][s,f]-norm),degree.ALPHA*a[1][s,f]))
        tails.append((degree.ALPHA*(b[0][s,f]-norm),degree.ALPHA*b[1][s,f]))
    paths,derivative,orientations=chain(regions,tails)
    means=derivative/degree.ALPHA
    source=degree.ALPHA*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
    output=-degree.ALPHA*pilot.LIMIT*math.log(degree.Z)
    correction=-degree.ALPHA*24575*math.log(degree.QUIET)
    score=(source+output+correction+math.log(4096)+paths)/LN2
    assert abs(score-11617.188069039783)<1e-6
    identity=float(means[0]+means[1]+means[2]-means[3])
    assert abs(identity-32767)<1e-6
    slopes={"restored_births_le_512":degree.ALPHA*(means[0]-512)/LN2,
            "quiet_le_24575":degree.ALPHA*(means[1]-24575)/LN2,
            "floor_4a_le_61503":degree.ALPHA*(4*means[2]-61503)/LN2}
    guard()
    print(json.dumps({"scope":"FLOAT positive-dual exact regional0/full16 submajorant; local branch sensitivity only",
        "whole_code_certificate":False,"full_source_class_bound":False,"directed":False,
        "parameters":{"z":degree.Z,"alpha":degree.ALPHA,"quiet_mark":degree.QUIET,"birth_mark":1,"floor_mark":1},
        "log2_control_score":score,"hierarchical_tilted_means":{"births_k":float(means[0]),
            "quiet_q":float(means[1]),"active_incoming_a":float(means[2]),"final_incoming16_probability":float(means[3])},
        "score_derivatives_bits_per_ln_raw_mark":{key:float(value) for key,value in slopes.items()},
        "counts_identity":{"q_plus_k_plus_a_minus_final":identity,"expected":32767},
        "interpretation":{"birth_negative_direction_helps_locally":bool(slopes["restored_births_le_512"]>0),
            "floor_negative_direction_helps_locally":bool(slopes["floor_4a_le_61503"]>0),
            "quiet_negative_direction_helps_locally":bool(slopes["quiet_le_24575"]>0),
            "joint_parameter_exhaustion_claim":False},
        "orientation_chain_bits":[value/LN2 for value,_ in orientations],"tests":{**checked,"route_normalizers":route_checks,
            "fixed_control_reproduction":"passed","full_count_identity":"passed"},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"exact_control_pin":PIN,
        "seconds":time.monotonic()-started,"peak_working_set_bytes":pool.peak_working_set(),
        "estimated_peak_bytes":128*(1<<20),"generated_bytes":0},indent=2))


if __name__=="__main__":
    main()
