"""Shared witnesses, outward binary64 adaptive updates, Arb final moments.

Every primitive in the nonnegative dynamic program rounds upward. Exact
power-of-two normalization keeps a shared exponent, never a floating log.
Arb constructs the initial region and Bernoulli coefficients. Search-side
hull pruning is not used: all bands are included in every maximum.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import general_occupancy as general
import tightened_occupancy as tight
import polynomial_regions as poly
import christoffel_caps as outer
from general_batch_certificate import density_cost


def up(value):
    return np.nextafter(value,np.inf)


def initial_float(region):
    # float(Arb) is nearest, followed by one successor. The Arb values are
    # singleton upper endpoints; nonzero underflow is rounded to min-subnormal.
    return up(np.array([[float(v) for v in row] for row in region])).reshape(-1,3,3)


def adaptive_matrices(region,left,right):
    current=region.copy();exponent=0
    for q in range(1,len(region)):
        a=current[:-1];b=current[1:]
        updated=up(up(left[0]*a)+up(right[0]*b))
        for x,y in zip(left[1:],right[1:]):
            np.maximum(updated,up(up(x*a)+up(y*b)),out=updated)
        maximum=float(updated.max());assert maximum>0 and math.isfinite(maximum)
        shift=math.frexp(maximum)[1]
        # Power-of-two rescaling is exact for normal results. nextafter also
        # covers any subnormal rounding in downward rescaling.
        current=up(np.ldexp(updated,-shift));exponent+=shift
        yield q,current[0].copy(),exponent


def build(screen_name,lower,upper,tag,verify=False):
    assert 1<=lower<=upper<=8192
    assert Path(screen_name).name==screen_name and tag.isidentifier()
    screen_path=base.HERE/'generated'/screen_name;output=base.HERE/'generated'/f'{tag}_outward.json'
    if not verify:assert not output.exists()
    old=base.read(output) if verify else None
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        for name,digest in old['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    screen=base.read(screen_path)
    for name,digest in screen['source_sha256'].items():assert base.sha(base.HERE/name)==digest
    bands=screen['bands'];assert sorted(w for band in bands for w in band)==list(base.WEIGHTS)
    outer.build(verify=True)
    caps=outer.deterministic_caps();t,s,spectrum=base.load_map(tight.NAME);kernel=tight.kernel_spectrum()
    ctx.prec=512 if verify else 256
    from audit_bch_q1_full_arb import rational
    best={};selected={}
    for index,witness in enumerate(screen['rows']):
        j=witness['witness_tenth'];assert isinstance(j,int) and -120<=j<=20
        lam=(arb(j)/10).exp();correction=(209716*lam).exp()
        region=initial_float(poly.regions(t,s,spectrum,kernel,(-lam).exp(),upper))
        ps=[base.decode(v) for v in witness['p']];assert len(ps)==len(bands) and all(0<p<1 for p in ps)
        costs=[density_cost(band,p,caps) for band,p in zip(bands,ps)]
        roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
        probabilities=[arb(p.numerator)/p.denominator for p in ps]
        left=up(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]))
        right=up(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]))
        for q,matrix,exponent in adaptive_matrices(region,left,right):
            if q<lower:continue
            value=tuple(arb(float(v)) for v in matrix.flat)
            for _ in range(8):value=q1.positive_mul(value,value)
            bound=math.comb(8192,q)*len(bands)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
            assert bound>0
            if q not in best or bound<best[q]:best[q]=bound;selected[q]=index
        passing=[q for q,u in best.items() if u<F(1,1<<60)]
        print('Shared witness',index,'anchor',witness['occupation'],'tilt',j,'occupancies below 2^-60',len(passing),'of',upper-lower+1,flush=True)
    total=sum(best.values(),F(0));rows=[]
    for q in range(lower,upper+1):
        bound=best[q]
        if old:assert bound<=base.decode(old['rows'][q-lower]['upper'])
        rows.append(dict(occupation=q,upper=base.encode(bound),witness_index=selected[q],
                         margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator)))
    if old:
        assert total<=base.decode(old['range_upper'])
        print('512-bit initial/final replay and full directed recurrence passed',flush=True);return
    local,outer_paths=general.q2.source_paths()
    local += [Path(__file__),Path(poly.__file__),Path(tight.__file__),Path(outer.__file__),Path(general.__file__),
              base.HERE/'general_batch_certificate.py',base.HERE/'certify_q3_compact.py',screen_path,
              base.HERE/'generated/christoffel_oa29_caps.json']
    local += [base.HERE/name for name in screen['source_sha256']]
    payload=dict(status='OUTWARD_SHARED_ADAPTIVE_RANGE',configuration=tight.NAME,occupancy_range=[lower,upper],
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,step_bits=t,state_bits=s,distance_cutoff=209716),
        bands=bands,rows=rows,range_upper=base.encode(total),range_below_2_to_minus_60=total<F(1,1<<60),
        all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer_paths})
    base.write_new(output,payload)
    print('Range',lower,upper,'margin',math.log2(total.denominator)-math.log2(total.numerator),
          'worst occupancies',sorted(((r['occupation'],round(r['margin_bits_diagnostic'],3)) for r in rows),key=lambda x:x[1])[:8],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',required=True);p.add_argument('--lower',type=int,required=True)
    p.add_argument('--upper',type=int,required=True);p.add_argument('--tag',required=True);p.add_argument('--verify',action='store_true')
    a=p.parse_args();build(a.screen,a.lower,a.upper,a.tag,a.verify)
