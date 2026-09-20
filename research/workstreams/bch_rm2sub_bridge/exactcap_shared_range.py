"""Full-occupancy shared bounds using the improved exact shell-cap snapshot.

Unlike the constant-row split, this includes the thirteenth, all-one band.
Each depth certifies all messages of that occupancy. Failed rows remain in
the receipt. The caller must check the range sum before claiming coverage.
"""
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
import screen_exponential_modes as cap_loader
import christoffel_caps as cap_base
from general_batch_certificate import density_cost


def run(lower,upper,tag,verify=False):
    assert 1<=lower<=upper<=8192 and tag.isidentifier()
    path=base.HERE/'generated'/f'exactcap_shared_{tag}_outward.json'
    screen_path=base.HERE/'generated/independent_p_anchors_screen.json'
    old=base.read(path) if verify else None
    if not verify:assert not path.exists()
    if old:
        assert old['occupancy_range']==[lower,upper]
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        caps={int(w):v for w,v in old['used_caps'].items()}
        sources=[base.HERE/name for name in old['cap_receipts']]
        rebuilt=cap_base.deterministic_caps()
        for source in sources:
            receipt=base.read(source)
            for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
            for name,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
            assert receipt['rational_primal_dual_checks_passed']
            w=receipt['weight'];rebuilt[w]=rebuilt[256-w]=min(rebuilt[w],receipt['cap'])
        assert caps==rebuilt
    else:caps,sources=cap_loader.latest_caps()
    screen=base.read(screen_path)
    for name,digest in screen['source_sha256'].items():assert base.sha(base.HERE/name)==digest
    bands=screen['bands'];assert sorted(w for band in bands for w in band)==list(base.WEIGHTS)
    witness=next(r for r in screen['rows'] if r['occupation']==1024)
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=512 if verify else 256
    lam=(arb(witness['witness_tenth'])/10).exp();correction=(209716*lam).exp()
    region=poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),upper)
    ps=[base.decode(v) for v in witness['p']]
    costs=[density_cost(band,p,caps) for band,p in zip(bands,ps)]
    roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probabilities=[arb(p.numerator)/p.denominator for p in ps]
    left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]),np.inf)
    right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]),np.inf)
    keep=hull.indices(left,right);print('Exact hull',len(keep),'of',len(bands),'bands',flush=True)
    mantissas,exponents=scaled.initial(region);rows=[];total=F(0)
    from audit_bch_q1_full_arb import rational
    for q,matrix,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep]):
        if q<lower:continue
        value=tuple(arb(float(v)) for v in matrix.flat)
        for _ in range(8):value=q1.positive_mul(value,value)
        bound=math.comb(8192,q)*len(bands)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
        assert bound>0;total+=bound
        if old:assert bound<=base.decode(old['rows'][q-lower]['upper'])
        rows.append(dict(occupation=q,upper=base.encode(bound),margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator)))
    if old:
        assert total<=base.decode(old['range_upper'])
        print('512-bit improved-cap shared range replay passed',flush=True);return
    local=[Path(__file__),Path(poly.__file__),Path(scaled.__file__),Path(hull.__file__),Path(tight.__file__),Path(q1.__file__),
        Path(cap_base.__file__),Path(cap_loader.__file__),base.HERE/'general_occupancy.py',base.HERE/'general_batch_certificate.py',
        base.HERE/'certify_q3_compact.py',screen_path]+sources+[base.HERE/name for name in screen['source_sha256']]
    base.write_new(path,dict(status='OUTWARD_EXACTCAP_SHARED_RANGE',configuration=tight.NAME,occupancy_range=[lower,upper],
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,step_bits=t,state_bits=s,distance_cutoff=209716),
        used_caps={str(w):caps[w] for w in base.WEIGHTS},cap_receipts=[str(p.relative_to(base.HERE)) for p in sources],
        rows=rows,range_upper=base.encode(total),range_below_2_to_minus_60=total<F(1,1<<60),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    passing=[r['occupation'] for r in rows if r['margin_bits_diagnostic']>60]
    print('Range margin',math.log2(total.denominator)-math.log2(total.numerator),'passing >60 bits',len(passing),
        'endpoints',(min(passing),max(passing)) if passing else None,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--lower',type=int,required=True);p.add_argument('--upper',type=int,required=True)
    p.add_argument('--tag',required=True);p.add_argument('--verify',action='store_true');a=p.parse_args()
    run(a.lower,a.upper,a.tag,a.verify)
