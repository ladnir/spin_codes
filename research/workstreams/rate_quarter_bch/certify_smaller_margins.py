"""Outward certification and higher-precision replay of quarter-rate margins.

Only fixed witness parameters come from binary64 discovery. All moments,
counting costs, powers, endpoint maxima, and unions are recomputed in Arb.
Recorded upper bounds are exact dyadics with upward-rounded 160-bit mantissas.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path

from flint import arb,arb_poly,ctx
import flint
import evaluate_fixed_inner as fixed
import smaller_outer
from verify_fixed_inner_results import check_coverage
import typed_dense_boxes as typed

B,L,N,T,S = 128,32768,4194304,128,19


def number(value):
    if isinstance(value,float): value=F.from_float(value)
    if isinstance(value,str): value=F(value)
    if isinstance(value,F): return arb(value.numerator)/value.denominator
    return arb(value)


def upper(value):
    assert value.is_finite(),'nonfinite interval'
    return value.upper()


def pack(value):
    point=upper(value)
    assert point>=0
    m,e=map(int,point.man_exp())
    shift=max(0,m.bit_length()-160)
    if shift: m=(m+(1<<shift)-1)>>shift; e+=shift
    return dict(mantissa=str(m),exponent=e)


def unpack(value):
    return arb(int(value['mantissa']))*arb(2)**int(value['exponent'])


def margin(value):
    return float(-value.log()/arb(2).log())


def identity():
    return tuple(arb(int(i==j)) for i in range(3) for j in range(3))


def mul(a,b):
    return tuple(a[3*i]*b[j]+a[3*i+1]*b[3+j]+a[3*i+2]*b[6+j] for i in range(3) for j in range(3))


def power(a,n):
    out=identity()
    while n:
        if n&1: out=mul(out,a)
        n>>=1
        if n: a=mul(a,a)
    return out


def terminal(a,n):
    return sum(power(a,n)[:3],arb(0))


def polynomial_mul(a,b,maximum):
    return tuple(sum((a[3*i+k]*b[3*k+j] for k in range(3)),arb_poly()).truncate(maximum+1)
                 for i in range(3) for j in range(3))


def epoch(spectrum,kernel,lam,maximum):
    den=(1<<S)-1; kappa=arb(den)/(den-1)
    z=(-lam).exp(); powers=[z**j for j in range(2*T+1)]
    rows=[]
    for j in range(maximum+1):
        total=math.comb(T,j)
        moments={w:upper(sum((math.comb(w,v)*math.comb(T-w,j-v)*powers[w+j-2*v]
                    for v in range(max(0,j-T+w),min(w,j)+1)),arb(0))/total) for w in spectrum}
        arbitrary=max(moments.values())
        uniform=upper(sum((n*moments[w] for w,n in spectrum.items()),arb(0))/den)
        nonkernel=arb(total-kernel[j])/total
        rows.append(tuple(map(upper,(arb(kernel[j])/total*powers[j],nonkernel*powers[j],arb(0),
                         min(upper(nonkernel),arbitrary)/den,arb(0),arbitrary,
                         kappa*min(upper(nonkernel),uniform)/den,arb(0),kappa*uniform))))
    return rows


def region(epoch_rows,maximum,linear=False):
    a=tuple(arb_poly([row[k]*math.comb(T,j) for j,row in enumerate(epoch_rows[:maximum+1])]) for k in range(9))
    out=tuple(arb_poly([int(i==j)]) for i in range(3) for j in range(3))
    n=L//T
    if linear:
        for _ in range(n): out=polynomial_mul(out,a,maximum)
    else:
        while n:
            if n&1: out=polynomial_mul(out,a,maximum)
            n>>=1
            if n: a=polynomial_mul(a,a,maximum)
    return [tuple(max(arb(0),upper(p[j]/math.comb(L,j))) for p in out) for j in range(maximum+1)]


def q1(spectrum,lam,linear=False):
    den=(1<<S)-1; kappa=arb(den)/(den-1); z=(-lam).exp(); d=min(spectrum)
    m0=sum((n*z**w for w,n in spectrum.items()),arb(0))/den
    m1=sum((n*(w*z**(w-1)+(T-w)*z**(w+1)) for w,n in spectrum.items()),arb(0))/(T*den)
    zero=(arb(1),arb(0),arb(0),arb(0),arb(0),z**d,arb(0),arb(0),kappa*m0)
    one=(arb(0),z,arb(0),z**(d-1)/den,arb(0),z**(d-1),kappa*m1/den,arb(0),kappa*m1)
    rz,ra=identity(),(arb(0),)*9
    if linear:
        for _ in range(L//T):
            ra=tuple(x+y for x,y in zip(mul(ra,zero),mul(rz,one))); rz=mul(rz,zero)
    else:
        n=L//T
        while n:
            if n&1:
                ra=tuple(x+y for x,y in zip(mul(ra,zero),mul(rz,one))); rz=mul(rz,zero)
            n>>=1
            if n:
                one=tuple(x+y for x,y in zip(mul(one,zero),mul(zero,one))); zero=mul(zero,zero)
    ra=tuple(v/(L//T) for v in ra)
    def vm(v,m):
        x,y,z=v
        return (x*m[0]+y*m[3]+z*m[6],x*m[1]+y*m[4]+z*m[7],x*m[2]+y*m[5]+z*m[8])
    current=[(arb(1),arb(0),arb(0))]
    for n in range(B):
        out=[(arb(0),)*3 for _ in range(n+2)]
        for w,v in enumerate(current):
            out[w]=tuple(x+y for x,y in zip(out[w],vm(v,rz)))
            out[w+1]=tuple(x+y for x,y in zip(out[w+1],vm(v,ra)))
        current=out
    return [sum(v,arb(0))/math.comb(B,w) for w,v in enumerate(current)]


def costs(counts,bands,ps):
    assert sorted(w for band in bands for w in band)==sorted(counts)
    result=[]
    for band,p in zip(bands,ps):
        if band==[B]:
            assert p==1
            result.append(arb(counts[B])); continue
        assert p>0 and p<1
        result.append(max(upper(arb(counts[w])/math.comb(B,w)/p**w/(1-p)**(B-w)) for w in band))
    return result


def parameters(counts,bands,shift):
    ps=[]
    for band in bands:
        if band==[B]: ps.append(arb(1)); continue
        midpoint=arb(min(band)+max(band))/(2*B)
        odds=midpoint/(1-midpoint)*number(shift).exp()
        ps.append(odds/(1+odds))
    return ps,costs(counts,bands,ps)


def mixture(rows,ps):
    probabilities=[arb(1)]
    for p in ps:
        out=[arb(0)]*(len(probabilities)+1)
        for j,v in enumerate(probabilities): out[j]+=v*(1-p); out[j+1]+=v*p
        probabilities=out
    return tuple(sum((p*rows[j][k] for j,p in enumerate(probabilities)),arb(0)) for k in range(9))


def adaptive(rows,ps,gammas,q):
    roots=[upper((g.log()/B).exp()) for g in gammas]
    current=rows[:q+1]
    for _ in range(q):
        current=[tuple(max(upper(root*((1-p)*left[k]+p*right[k])) for root,p in zip(roots,ps))
                       for k in range(9)) for left,right in zip(current[:-1],current[1:])]
    return terminal(current[0],B)*math.comb(L,q)*len(ps)**q


def dense_box(box,banks,bands,counts,epoch_at,cutoff):
    witness=box['witness']; raw=banks[witness['probability_bank']]
    ps=[number(v) for v in raw['probabilities']]
    assert ps[0]==0
    gammas=[arb(1)]+costs(counts,bands,ps[1:])
    fractions=[F.from_float(v) for v in witness['proposal']]
    assert all(v>0 for v in fractions)
    total=sum(fractions); proposal=[number(v/total) for v in fractions]
    theta=sum((p*v for p,v in zip(ps,proposal)),arb(0))
    z=F.from_float(witness['log_surprisal']); lam=number(z).exp()
    rows=epoch_at(z)
    probabilities=[math.comb(T,j)*theta**j*(1-theta)**(T-j) for j in range(T+1)]
    matrix=tuple(sum((p*rows[j][k] for j,p in enumerate(probabilities)),arb(0)) for k in range(9))
    moment=terminal(matrix,N//T).log()
    corners=typed.vertices(box['lower'],box['upper'],L)
    values=[]
    for corner in corners:
        c=list(map(int,corner)); mult=math.comb(L,c[0])*math.comb(L-c[0],c[1])
        value=cutoff*lam+moment-(B-1)*arb(mult).log()
        value+=sum((n*(g.log()-B*p.log()) for n,g,p in zip(c,gammas,proposal)),arb(0))
        values.append(upper(value))
    widths=[hi-lo+1 for lo,hi in zip(box['lower'],box['upper'])]
    lattice=math.prod(widths)//max(widths)
    return max(values).exp()*lattice


def run(precision,verify):
    ctx.prec=precision
    witness_path=fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json'
    receipt_path=fixed.HERE/'SMALLER_MARGIN_CERTIFICATE.json'
    witnesses=json.loads(witness_path.read_text())
    for name,digest in witnesses['source_sha256'].items():
        assert fixed.sha(fixed.ROOT/name)==digest,name
    saved=json.loads(receipt_path.read_text()) if verify else None
    if saved:
        for name,digest in saved['source_sha256'].items(): assert fixed.sha(fixed.ROOT/name)==digest,name
        assert precision>saved['precision_bits']
    spectrum,kernel,sources=fixed.load_inner()
    construction=smaller_outer.construction()
    counts={w:n for w,n in enumerate(smaller_outer.spectrum()) if w and n}
    assert witnesses['outer']==[B,32,32] and witnesses['message_exponent']==20
    assert [r['distance_target'] for r in witnesses['results']]==['33/200','19/100']
    if saved:
        assert saved['status']=='OUTWARD_MARGIN_CERTIFICATE'
        assert saved['outer']==construction and saved['message_bits']==1<<20 and saved['output_bits']==N
        assert saved['inner']==dict(t=T,s=S) and len(saved['results'])==len(witnesses['results'])
    @lru_cache(maxsize=256)
    def epoch_at(z): return epoch(spectrum,kernel,number(z).exp(),T)
    @lru_cache(maxsize=256)
    def region_at(z,maximum): return region(epoch_at(z),maximum,linear=verify and maximum<=2)
    @lru_cache(maxsize=256)
    def q1_at(z): return q1(spectrum,number(z).exp(),linear=verify)
    @lru_cache(maxsize=32)
    def params(kind,shift): return parameters(counts,witnesses[kind+'_bands'],shift)
    results=[]
    for index,row in enumerate(witnesses['results']):
        old=saved['results'][index] if saved else None
        cutoff=row['bad_weight']; delta=F(row['distance_target']); target=40 if delta==F(33,200) else 30
        assert cutoff==N*delta.numerator//delta.denominator
        if old:
            assert old['distance_target']==row['distance_target'] and old['bad_weight']==cutoff
            assert old['target_bits']==target and old['maximum_sparse']==row['maximum_sparse']
            assert old['exact_occupation_coverage_checked'] is True
        def record(value,prior=None):
            if prior is not None:
                assert upper(value)<=unpack(prior),'higher-precision bound exceeds certificate'
                return prior
            return pack(value)
        q1_terms=[]
        assert sorted(w['weight'] for w in row['q1'])==sorted(counts)
        for j,w in enumerate(row['q1']):
            z=F(w['log_tilt']); lam=number(z).exp()
            value=L*counts[w['weight']]*q1_at(z)[w['weight']]*(cutoff*lam).exp()
            q1_terms.append(record(value,old['q1_terms'][j] if old else None))
        pair_terms=[]
        expected=[list(map(int,c)) for c in fixed.composition.compositions(len(witnesses['q2_bands']),2)]
        assert [w['band_indices'] for w in row['q2']]==expected
        for j,w in enumerate(row['q2']):
            z=F(w['log_tilt']); ps,gs=params('q2',w['shift']); ids=w['band_indices']
            mult=2 if ids[0]!=ids[1] else 1
            mat=mixture(region_at(z,2),[ps[g] for g in ids])
            value=terminal(mat,B)*gs[ids[0]]*gs[ids[1]]*mult*math.comb(L,2)*(cutoff*number(z).exp()).exp()
            pair_terms.append(record(value,old['q2_terms'][j] if old else None))
        adaptive_terms=[]
        assert [w['occupation'] for w in row['adaptive']]==list(range(3,row['maximum_sparse']+1))
        for j,w in enumerate(row['adaptive']):
            z=F(w['log_tilt']); q=w['occupation']; ps,gs=params('adaptive',w['shift'])
            value=adaptive(region_at(z,q),ps,gs,q)*(cutoff*number(z).exp()).exp()
            adaptive_terms.append(record(value,old['adaptive_terms'][j] if old else None))
            if q%16==0: print(precision,row['distance_target'],'sparse',q,flush=True)
        dense=row['dense']; minimum=row['maximum_sparse']+1
        assert dense['occupation_min']==minimum and dense['occupation_max']==L
        check_coverage(dense['selected_boxes'],L,minimum)
        dense_terms=[]
        for j,box in enumerate(dense['selected_boxes']):
            value=dense_box(box,dense['probability_banks'],dense['bands'],counts,epoch_at,cutoff)
            dense_terms.append(record(value,old['dense_terms'][j] if old else None))
            if j%32==0: print(precision,row['distance_target'],'dense',j+1,'/',len(dense['selected_boxes']),flush=True)
        groups=dict(q1_terms=q1_terms,q2_terms=pair_terms,adaptive_terms=adaptive_terms,dense_terms=dense_terms)
        bounds={name:record(sum((unpack(v) for v in terms),arb(0)),old['group_upper'][name] if old else None)
                for name,terms in groups.items()}
        total=record(sum((unpack(v) for v in bounds.values()),arb(0)),old['union_upper'] if old else None)
        assert unpack(total)<arb(2)**(-target),'target not certified'
        result=dict(distance_target=row['distance_target'],bad_weight=cutoff,target_bits=target,
                    maximum_sparse=row['maximum_sparse'],**groups,group_upper=bounds,union_upper=total,
                    margin_bits_diagnostic=margin(unpack(total)),exact_occupation_coverage_checked=True)
        results.append(result)
        print('PASS',precision,row['distance_target'],result['margin_bits_diagnostic'],flush=True)
    sources += [Path(__file__),witness_path,Path(fixed.__file__),Path(smaller_outer.__file__),Path(fixed.outer.__file__),
                Path(fixed.composition.__file__),Path(typed.__file__),fixed.HERE/'verify_fixed_inner_results.py',
                fixed.HERE/'sources/EBCH128_29.wd',fixed.HERE/'sources/EBCH128_36.wd']
    if verify:
        output=fixed.HERE/'SMALLER_MARGIN_CERTIFICATE_REPLAY.json'
        payload=dict(status='HIGHER_PRECISION_OUTWARD_REPLAY_PASSED',precision_bits=precision,
                     producer_sha256=fixed.sha(receipt_path),source_sha256=saved['source_sha256'],
                     margins_bits=[r['margin_bits_diagnostic'] for r in results],
                     q1_and_q2_region_replay='linear epoch iteration; producer uses binary powering')
    else:
        output=receipt_path
        payload=dict(status='OUTWARD_MARGIN_CERTIFICATE',precision_bits=precision,flint_version=flint.__version__,
                     outer=construction,message_bits=1<<20,output_bits=N,inner=dict(t=T,s=S),
                     results=results,source_sha256={p.relative_to(fixed.ROOT).as_posix():fixed.sha(p) for p in sources})
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--precision',type=int)
    args=parser.parse_args()
    run(args.precision or (512 if args.verify else 256),args.verify)
