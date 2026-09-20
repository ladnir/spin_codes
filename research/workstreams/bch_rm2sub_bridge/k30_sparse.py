"""Refresh-aware sparse occupancies at K=2^30, using exact kernel exclusion.

For Q<=7 every nonempty epoch input has nonzero syndrome (kernel distance 8).
No larger-Q claim is made by this module.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, arb_poly, ctx
import bridge as base
import larger_state_maps as maps
import certify_refresh_q1 as matrix
from certify_k30_q1 import ROWS,CUTOFF
from occupation_three import BANDS
from audit_bch_q1_full_arb import rational
from general_batch_certificate import density_cost
from migrate_legacy_workspace import restore_or_verify


def caps():
    data=base.read(base.HERE/'generated/christoffel_oa29_caps.json')
    result={row['weight']:row['cap'] for row in data['rows']}
    for path in sorted((base.HERE/'generated').glob('joint_shell_*/cap.json')):
        row=base.read(path);w=row['weight']
        result[w]=result[256-w]=min(result[w],row['cap'])
    return result


def epochs(t,s,spectrum,kernel,z,maximum,number):
    assert 0<=maximum<min(w for w,c in kernel.items() if w and c)
    result=[matrix.epoch(t,s,spectrum,z,number)[0]]
    m=(1<<s)-1;kappa=number(m)/(m-1)
    for j in range(1,maximum+1):
        uniform=sum((number(n)*sum((number(math.comb(w,v)*math.comb(t-w,j-v))*z**(w+j-2*v)
                    for v in range(max(0,j-t+w),min(w,j)+1)),number(0))/math.comb(t,j)
                    for w,n in spectrum.items()),number(0))/m
        arbitrary=z**max(0,min(abs(w-j) for w in spectrum))
        row=[number(0)]*16;row[1]=z**j
        for state,moment in ((1,arbitrary),(2,uniform),(3,kappa*uniform)):
            row[4*state]=moment/m;row[4*state+3]=moment*(m-1)/m
        result.append(tuple(row))
    return result


def regions(t,s,spectrum,kernel,z,maximum,length=ROWS,reverse=False):
    assert length%t==0
    epoch=epochs(t,s,spectrum,kernel,z,maximum,arb)
    power=tuple(arb_poly([row[k]*math.comb(t,j) for j,row in enumerate(epoch)]) for k in range(16))
    current=tuple(arb_poly([int(i==j)]) for i in range(4) for j in range(4))
    def mul(a,b):
        return tuple(sum((a[4*i+k]*b[4*k+j] for k in range(4)),arb_poly()).truncate(maximum+1)
                     for i in range(4) for j in range(4))
    count=length//t
    if reverse:
        for bit in bin(count)[2:]:
            current=mul(current,current)
            if bit=='1':current=mul(current,power)
    else:
        while count:
            if count&1:current=mul(current,power)
            count>>=1
            if count:power=mul(power,power)
    return [tuple(max(arb(0),(p[j]/math.comb(length,j)).upper()) for p in current)
            for j in range(maximum+1)]


def terminal_float(m):
    scale=0.
    for _ in range(8):
        m=m@m;factor=float(m.max());m/=factor;scale=2*scale+math.log(factor)
    return math.log(float(m[0].sum()))+scale


def choose_ps(region,q,counts):
    ps=[]
    for band in BANDS[:-1]:
        def objective(theta):
            p=1/(1+math.exp(-theta))
            cost=max(math.log(counts[w])-math.log(math.comb(256,w))-w*math.log(p)-(256-w)*math.log1p(-p) for w in band)
            weights=np.array([math.comb(q,j)*p**j*(1-p)**(q-j) for j in range(q+1)])
            return terminal_float(np.einsum('i,ijk->jk',weights,region))+q*cost
        fit=minimize_scalar(objective,bounds=(-7,8),method='bounded')
        assert fit.success
        ps.append(F.from_float(1/(1+math.exp(-float(fit.x)))))
    return ps+[F(1)]


def adaptive(region,q,ps,counts,number):
    roots=[]
    for band,p in zip(BANDS,ps):
        cost=F(1) if p==1 else density_cost(band,p,counts)
        if number is float:roots.append(math.exp(math.log(cost)/256))
        else:roots.append(((arb(cost.numerator)/cost.denominator).log()/256).exp().upper())
    if number is float:
        current=np.array(region,dtype=float).reshape(q+1,4,4)
        for _ in range(q):
            current=np.maximum.reduce([root*((1-float(p))*current[:-1]+float(p)*current[1:]) for root,p in zip(roots,ps)])
        return current[0]
    probabilities=[arb(p.numerator)/p.denominator for p in ps]
    current=region
    for _ in range(q):
        current=[tuple(max((root*((1-p)*current[j][k]+p*current[j+1][k])).upper()
                           for root,p in zip(roots,probabilities)) for k in range(16))
                 for j in range(len(current)-1)]
    return current[0]


def run(output,verify=False):
    output=output.resolve();saved=base.read(output) if verify else None
    if not verify:assert not output.exists()
    restore_or_verify(base.ROOT)
    t,s,spectrum,kernel=maps.load('t64_s20');counts=caps()
    ctx.prec=512 if verify else 256
    results=[]
    for q in range(2,8):
        if verify:
            old=next(row for row in saved['rows'] if row['occupation']==q)
            candidates=[(base.decode(old['scaled_tilt']),[base.decode(p) for p in old['p']])]
        else:
            candidates=[(F(q*n,2),None) for n in (3,5,8,12,18,26,40)]
        best=None
        for a,ps in candidates:
            lam=(arb(a.numerator)/a.denominator)/ROWS
            region=regions(t,s,spectrum,kernel,(-lam).exp(),q,reverse=verify)
            if ps is None:
                approx=np.array([[float(v) for v in m] for m in region]).reshape(q+1,4,4)
                ps=choose_ps(approx,q,counts)
            m=adaptive(region,q,ps,counts,arb)
            for _ in range(8):m=matrix.product(m,m)
            upper=math.comb(ROWS,q)*len(BANDS)**q*rational((sum(m[:4],arb(0))*(CUTOFF*lam).exp()).upper())
            if best is None or upper<best[0]:best=(upper,a,ps)
        bound,a,ps=best
        assert bound<F(1,1<<44),f'Q{q} misses budget'
        if verify:
            assert bound<=base.decode(old['upper'])
        results.append(dict(occupation=q,upper=base.encode(bound),scaled_tilt=base.encode(a),
            p=[base.encode(p) for p in ps],margin_bits=math.log2(bound.denominator)-math.log2(bound.numerator)))
        print('replay' if verify else 'producer','Q',q,'margin',results[-1]['margin_bits'],flush=True)
    if verify:
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        assert sum((base.decode(r['upper']) for r in saved['rows']),F(0))==base.decode(saved['range_upper'])
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='K30_Q2_THROUGH_Q7_512_BIT_REPLAY_PASSED',producer_sha256=base.sha(output),full_distance_proved=False))
    else:
        paths=[Path(__file__),Path(matrix.__file__),base.HERE/'certify_k30_q1.py',
               base.HERE/'general_batch_certificate.py',base.HERE/'occupation_three.py',base.HERE/'MIGRATION_MANIFEST.json']
        base.write_new(output,dict(status='K30_SPARSE_OUTWARD_CERTIFICATE',
            parameters=dict(message_bits=1<<30,outer_rows=ROWS,cutoff=CUTOFF,step_bits=64,state_bits=20),
            rows=results,range_upper=base.encode(sum((base.decode(r['upper']) for r in results),F(0))),
            full_distance_proved=False,source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();run(args.output,args.verify)
