"""Derive a dimension-69 hull anchor from the published dimension-71 spectrum."""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from audit_bch_hulls import rows,hull,contains
from audit_bch_quadratic_sums import signed_sum
from audit_bch_closure_envelope import kraw_table
from anchors import L71_HALF
from prepare_bch_hull_probe import build as previous_build,export
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_anchor_probe'


def build():
    l=rows(59,71)
    r=hull(l)
    hq=hull(rows(39,123))
    assert len(r)==69 and all(contains(hq,x) for x in r)
    gauss=signed_sum(l)
    assert gauss['bilinear_rank']==2 and int(gauss['signed_sum'])==-(1<<70)
    assert all(x.bit_count()%4==0 for x in r)
    # Exactly 2^69 words of L have weight divisible by four. Its doubly-even
    # hull already has that size, so those words are precisely its hull.
    anchor=[L71_HALF.get(min(w,256-w),0) if w%4==0 else 0 for w in range(257)]
    assert sum(anchor)==1<<69
    kt=kraw_table()
    dual=[]
    for row in kt:
        value,rem=divmod(sum(x*y for x,y in zip(anchor,row)),1<<69)
        assert rem==0 and value>=0
        dual.append(value)
    assert sum(dual)==1<<187
    model,scales=previous_build()
    for j in range(0,129,2):
        if j%4==0:
            model['constraints'].append(dict(name='HQ_ge_Lhull_'+str(j),
                coeffs={f'r_{j}':'1',f's_{j}':'255'},sense='ge',rhs=str(anchor[j])))
        coeffs={}
        for a,mult in (('r',1),('s',255)):
            for w in range(0,129,4):
                coeffs[f'{a}_{w}']=str(mult*(kt[j][w]+(kt[j][256-w] if w!=128 else 0)))
        model['constraints'].append(dict(name='HQdual_le_Lhulldual_'+str(j),
                coeffs=coeffs,sense='le',rhs=str((1<<93)*dual[j])))
    model['metadata']['Lhull_anchor']={
        'dimension':69,'hull_basis_hex':[hex(x) for x in r],
        'half_spectrum':{str(w):str(anchor[w]) for w in range(0,129,4)},
        'dual_half_spectrum':{str(w):str(dual[w]) for w in range(0,129,2)},
        'gauss':gauss,'exact_HQ_containment_verified':True,
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/anchors.py',ROOT/'code/audit_bch_hulls.py',ROOT/'code/prepare_bch_hull_probe.py')}}
    return model,scales


if __name__=='__main__':
    model,scales=build()
    lp,count=export(model,scales)
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'scales.json',scales)
    with (FOLDER/'h_38.lp').open('x') as f:
        f.write(lp)
    print(json.dumps(dict(variables=len(scales),rows=count)))
