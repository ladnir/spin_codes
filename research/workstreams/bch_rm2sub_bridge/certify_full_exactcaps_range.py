"""Outward shared witnesses for complete occupancies with captured exact caps."""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import tightened_occupancy as tight
import polynomial_regions as poly
import scaled_adaptive as scaled
import positive_line_hull as hull
import christoffel_caps as cap_base
from general_batch_certificate import density_cost


def checked_caps(screen):
    caps=cap_base.deterministic_caps()
    for name,digest in screen['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for name in screen['cap_receipts']:
        receipt=base.read(base.HERE/name)
        for p,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/p)==digest
        for p,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/p)==digest
        assert receipt['rational_primal_dual_checks_passed']
        w=receipt['weight'];caps[w]=caps[256-w]=min(caps[w],receipt['cap'])
    assert set(map(int,screen['used_caps']))==set(base.WEIGHTS)
    assert all(caps[w]==screen['used_caps'][str(w)] for w in base.WEIGHTS)
    return caps


def run(screen_name,lower,upper,tag,verify=False):
    assert 1<=lower<=upper<=8192 and tag.isidentifier() and Path(screen_name).name==screen_name
    path=base.HERE/'generated'/f'full_exactcaps_{tag}_outward.json';screen_path=base.HERE/'generated'/screen_name
    old=base.read(path) if verify else None
    if not verify:assert not path.exists()
    if old:
        assert old['occupancy_range']==[lower,upper]
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    screen=base.read(screen_path);caps=checked_caps(screen);bands=screen['bands']
    assert sorted(w for band in bands for w in band)==list(base.WEIGHTS)
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=512 if verify else 256
    best={};chosen={};cached={}
    from audit_bch_q1_full_arb import rational
    for index,witness in enumerate(screen['rows']):
        tilt=witness['witness_tenth']
        if tilt not in cached:
            lam=(arb(tilt)/10).exp()
            cached[tilt]=(poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),upper),(209716*lam).exp())
        region,correction=cached[tilt]
        ps=[base.decode(v) for v in witness['p']];assert len(ps)==len(bands) and all(0<p<1 for p in ps)
        costs=[density_cost(band,p,caps) for band,p in zip(bands,ps)]
        roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
        probabilities=[arb(p.numerator)/p.denominator for p in ps]
        left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]),np.inf)
        right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]),np.inf)
        keep=hull.indices(left,right);mantissas,exponents=scaled.initial(region)
        for q,matrix,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep]):
            if q<lower:continue
            value=tuple(arb(float(v)) for v in matrix.flat)
            for _ in range(8):value=q1.positive_mul(value,value)
            bound=math.comb(8192,q)*len(bands)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
            assert bound>0
            if q not in best or bound<best[q]:best[q]=bound;chosen[q]=index
        print('Witness',index,'anchor',witness['occupation'],'tilt',tilt,'cumulative passing',sum(v<F(1,1<<60) for v in best.values()),flush=True)
    rows=[dict(occupation=q,upper=base.encode(best[q]),witness_index=chosen[q],
        margin_bits_diagnostic=math.log2(best[q].denominator)-math.log2(best[q].numerator)) for q in range(lower,upper+1)]
    total=sum(best.values(),F(0))
    if old:
        assert len(rows)==len(old['rows'])
        for row,saved in zip(rows,old['rows']):
            assert row['occupation']==saved['occupation'] and base.decode(row['upper'])<=base.decode(saved['upper'])
        assert total<=base.decode(old['range_upper'])
        print('512-bit full shared range replay passed',flush=True);return
    local=[Path(__file__),Path(poly.__file__),Path(scaled.__file__),Path(hull.__file__),Path(tight.__file__),Path(q1.__file__),
        Path(cap_base.__file__),base.HERE/'general_occupancy.py',base.HERE/'general_batch_certificate.py',
        base.HERE/'certify_q3_compact.py',screen_path]+[base.HERE/name for name in screen['local_sha256']]
    base.write_new(path,dict(status='OUTWARD_FULL_EXACTCAPS_RANGE',configuration=tight.NAME,occupancy_range=[lower,upper],
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,step_bits=t,state_bits=s,distance_cutoff=209716),
        rows=rows,range_upper=base.encode(total),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    passing=[r['occupation'] for r in rows if r['margin_bits_diagnostic']>60]
    print('Range margin',math.log2(total.denominator)-math.log2(total.numerator),'passing',len(passing),
        'endpoints',(min(passing),max(passing)) if passing else None,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',required=True);p.add_argument('--lower',type=int,required=True)
    p.add_argument('--upper',type=int,required=True);p.add_argument('--tag',required=True);p.add_argument('--verify',action='store_true');a=p.parse_args()
    run(a.screen,a.lower,a.upper,a.tag,a.verify)
