"""All-weight RM2Sub kernel envelope and general-occupancy band majorant.

The formulas apply to every 1<=Q<=8192. Computations truncate region
polynomials at the requested Q; numerical coverage is recorded separately.
"""
import argparse
import itertools
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import occupation_two as q2
from certify_q3_compact import dyadic_upper

NAME='t128_s15'
BANDS=(tuple(range(38,56,2)),tuple(range(56,88,2)),tuple(range(88,130,2)),
       tuple(range(130,186,2)),tuple(range(186,220,2))+(256,))


def kernel_spectrum():
    rows=base.read(base.HERE/'inputs'/f'{NAME}_b_kernel_spectrum.json')['by_total_weight']
    return {row['total_weight']:row['kernel_words'] for row in rows}


def epoch_matrices(t,s,spectrum,kernel,z,number,maximum):
    den=(1<<s)-1; kappa=number(den)/(den-1); out=[]
    for j in range(min(maximum,t)+1):
        beta=number(kernel.get(j,0))/math.comb(t,j)
        moments=[]
        for w,c in spectrum.items():
            value=sum((number(math.comb(w,v)*math.comb(t-w,j-v))*z**(w+j-2*v)
                       for v in range(max(0,j-t+w),min(w,j)+1)),number(0))/math.comb(t,j)
            moments.append((value,c))
        uniform=sum((value*c for value,c in moments),number(0))/den
        # Maximum over actual A weights bounds every arbitrary nonzero state.
        # Arb interval maxima are not ordered: an explicit upper endpoint is used.
        arbitrary=max(float(v) for v,c in moments) if number is float else None
        if number is F: arbitrary=max(v for v,c in moments)
        if number is arb:
            from audit_bch_q1_full_arb import rational
            upper=max(rational(v.upper()) for v,c in moments)
            arbitrary=arb(upper.numerator)/upper.denominator
        nonkernel=1-beta
        # Weightwise pointwise moment <=1. For termination, either marginal
        # upper bound can be used; taking the sum-free choice nonkernel is safe.
        terminate_d=nonkernel/den
        terminate_l=kappa*nonkernel/den
        if j==0:terminate_d=terminate_l=number(0)
        # A second valid bound is the full emission moment divided by M.
        # Use deterministic kernel counts to choose the zero-syndrome cases;
        # otherwise the moment bound avoids interval min operations.
        if kernel.get(j,0)==0:
            terminate_d=arbitrary/den;terminate_l=kappa*uniform/den
        out.append((beta*z**j,nonkernel*z**j,number(0),
                    terminate_d,number(0),arbitrary,
                    terminate_l,number(0),kappa*uniform))
    return out


def regions(t,s,spectrum,kernel,z,number,maximum,length=8192):
    assert 0<=maximum<=length and length%t==0
    epoch=epoch_matrices(t,s,spectrum,kernel,z,number,maximum)
    weighted=[tuple(v*math.comb(t,j) for v in m) for j,m in enumerate(epoch)]
    current=[tuple(number(int(i==j)) for i in range(3) for j in range(3))]
    for step in range(length//t):
        limit=min(maximum,(step+1)*t);updated=[]
        for degree in range(limit+1):
            terms=[q1.positive_mul(current[degree-j],weighted[j])
                   for j in range(max(0,degree-len(current)+1),min(degree,len(weighted)-1)+1)]
            updated.append(tuple(sum((v[k] for v in terms),number(0)) for k in range(9)))
        current=updated
    return [tuple(v/math.comb(length,j) for v in current[j]) for j in range(maximum+1)]


def compositions(q,k=5):
    if k==1:
        yield (q,);return
    for count in range(q+1):
        for rest in compositions(q-count,k-1):yield (count,)+rest


def multiplicity(counts):
    value=math.factorial(sum(counts))
    for count in counts:value//=math.factorial(count)
    return value


def distribution(counts,ps,number):
    current=[number(1)]
    for count,p in zip(counts,ps):
        p=number(p.numerator)/p.denominator if isinstance(p,F) else number(p)
        for _ in range(count):
            updated=[number(0)]*(len(current)+1)
            for j,v in enumerate(current):updated[j]+=v*(1-p);updated[j+1]+=v*p
            current=updated
    return current


def log_power(matrix):
    scale=0.
    for _ in range(8):
        matrix=matrix@matrix
        factor=float(matrix.max());assert factor>0
        matrix/=factor;scale=2*scale+math.log(factor)
    return math.log(float(matrix[0].sum()))+scale


def screen(occupancies):
    t,s,spectrum=base.load_map(NAME);kernel=kernel_spectrum();caps=q2.deterministic_caps()
    assert sorted(w for band in BANDS for w in band)==list(base.WEIGHTS)
    weights=[np.array(band) for band in BANDS]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in BANDS]
    midpoint=[(band[0]+band[-1])/512 for band in BANDS]
    selected={q:{c:(math.inf,None,None) for c in compositions(q)} for q in occupancies}
    for tenth in (-80,-75,-70,-65,-60,-55):
        lam=math.exp(tenth/10)
        region=np.array(regions(t,s,spectrum,kernel,math.exp(-lam),float,max(occupancies))).reshape(-1,3,3)
        for q in occupancies:
            for counts in selected[q]:
                def objective(shift,ret=False):
                    factor=math.exp(shift)
                    ps=[p*factor/(1-p+p*factor) for p in midpoint]
                    cost=sum(count*float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                             for count,v,w,p in zip(counts,logs,weights,ps) if count)
                    dist=np.asarray(distribution(counts,ps,float))
                    moment=log_power(np.sum(dist[:,None,None]*region[:q+1],axis=0))
                    value=moment+209716*lam+cost
                    return (value,ps) if ret else value
                opt=minimize_scalar(objective,bounds=(-2.,8.),method='bounded',options={'xatol':0.003})
                value,ps=objective(float(opt.x),True)
                if value<selected[q][counts][0]:selected[q][counts]=(value,tenth,ps)
        print('General occupancy tilt',tenth,'complete',flush=True)
    rows=[]
    for q,best in selected.items():
        log=base.lse([v[0]+math.log(multiplicity(c)) for c,v in best.items()])+math.log(math.comb(8192,q))
        rows.append(dict(occupation=q,margin_bits_diagnostic=-log/math.log(2),
            compositions=[dict(counts=c,witness_tenth=v[1],p=[base.encode(F.from_float(p)) for p in v[2]]) for c,v in best.items()]))
        print('Q',q,'screen margin',-log/math.log(2),flush=True)
    path=base.HERE/'generated'/f'general_q{min(occupancies)}_q{max(occupancies)}_screen.json'
    base.write_new(path,dict(status='GENERAL_OCCUPANCY_SCREEN_ONLY',bands=BANDS,rows=rows,
        source_sha256={'general_occupancy.py':base.sha(Path(__file__))}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--occupancies',type=int,nargs='+',default=[4,5,6,7,8])
    args=parser.parse_args();assert all(1<=q<=8192 for q in args.occupancies)
    screen(args.occupancies)
