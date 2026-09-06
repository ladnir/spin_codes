"""Outward certificates for specified (ordinary rows, all-one rows) cells.

This does not certify the unlisted cells of the two-dimensional range.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import tightened_occupancy as tight
import polynomial_regions as poly
import scaled_adaptive as scaled
from general_batch_certificate import density_cost


def run(screen_name,tag,verify=False):
    assert Path(screen_name).name==screen_name and tag.isidentifier()
    path=base.HERE/'generated'/screen_name;output=base.HERE/'generated'/f'constant_split_{tag}_outward.json'
    old=base.read(output) if verify else None
    if not verify:assert not output.exists()
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    screen=base.read(path)
    for name,digest in screen['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    bands=screen['bands'];assert sorted(w for band in bands for w in band)==list(base.WEIGHTS[:-1])
    caps={int(w):v for w,v in screen['used_caps'].items()}
    # Reconstruct the selected minimum from immutable cap receipts. Later
    # optimizer outputs cannot silently change an already selected witness.
    import christoffel_caps as cap_base
    reconstructed=cap_base.deterministic_caps()
    for name in screen['local_sha256']:
        if name.endswith('cap.json'):
            receipt=base.read(base.HERE/name)
            for relative,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/relative)==digest
            for relative,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/relative)==digest
            assert receipt['rational_primal_dual_checks_passed']
            w=receipt['weight'];reconstructed[w]=reconstructed[256-w]=min(reconstructed[w],receipt['cap'])
    assert all(caps[w]==reconstructed[w] for w in base.WEIGHTS)
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=512 if verify else 256;cached={};rows=[];total=F(0)
    from audit_bch_q1_full_arb import rational
    for index,row in enumerate(screen['rows']):
        d=row['ordinary_rows'];h=row['all_one_rows'];tilt=row['witness_tenth']
        assert 0<=d and 0<=h and 1<=d+h<=8192
        if tilt not in cached:
            lam=(arb(tilt)/10).exp()
            cached[tilt]=(poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),8192),(209716*lam).exp())
        region,correction=cached[tilt]
        if d:
            ps=[base.decode(v) for v in row['p']];assert len(ps)==len(bands) and all(0<p<1 for p in ps)
            gamma=[density_cost(band,p,caps) for band,p in zip(bands,ps)]
            roots=[((arb(v.numerator)/v.denominator).log()/256).exp().upper() for v in gamma]
            probabilities=[arb(p.numerator)/p.denominator for p in ps]
            left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]),np.inf)
            right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]),np.inf)
            mantissas,exponents=scaled.initial(region[h:h+d+1])
            for depth,matrix,exponent in scaled.matrices(mantissas,exponents,left,right):
                if depth%2048==0:print('Directed cell',d,h,'depth',depth,flush=True)
            assert depth==d
            value=tuple(arb(float(v)) for v in matrix.flat)
        else:value=region[h];exponent=0
        for _ in range(8):value=q1.positive_mul(value,value)
        bound=(math.comb(8192,d)*math.comb(8192-d,h)*len(bands)**d
               *rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper()))
        assert bound>0
        if old:assert bound<=base.decode(old['rows'][index]['upper'])
        total+=bound
        rows.append(dict(ordinary_rows=d,all_one_rows=h,occupation=d+h,upper=base.encode(bound),
                         margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator)))
        print('Outward constant split',d,h,'margin',rows[-1]['margin_bits_diagnostic'],flush=True)
    if old:
        assert total<=base.decode(old['listed_cells_upper'])
        print('512-bit input/final replay and directed cell recurrences passed',flush=True);return
    local=[Path(__file__),Path(poly.__file__),Path(tight.__file__),Path(scaled.__file__),Path(q1.__file__),Path(cap_base.__file__),
           base.HERE/'general_occupancy.py',base.HERE/'general_batch_certificate.py',base.HERE/'certify_q3_compact.py',path]
    local += [base.HERE/name for name in screen['local_sha256']]
    base.write_new(output,dict(status='OUTWARD_LISTED_CONSTANT_SPLIT_CELLS',configuration=tight.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,step_bits=t,state_bits=s,distance_cutoff=209716),
        rows=rows,listed_cells_upper=base.encode(total),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',required=True);p.add_argument('--tag',required=True);p.add_argument('--verify',action='store_true')
    a=p.parse_args();run(a.screen,a.tag,a.verify)
