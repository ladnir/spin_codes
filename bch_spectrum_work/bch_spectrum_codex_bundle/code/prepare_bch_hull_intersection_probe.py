"""Couple J=HP intersect hull(L71) to both spectra and their duals."""
import json
import math
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from audit_bch_hulls import rows,hull,nullspace,echelon,contains,cyclic_description
from bch_quotient import binary_poly_divmod,multiply_by_x_mod
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_anchor_probe import build as previous_build
from prepare_bch_hull_probe import export
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_intersection_probe'


def build():
    hp=hull(rows(37,131))
    hq=hull(rows(39,123))
    anchor=hull(rows(59,71))
    # (A intersect B)^perp = A^perp + B^perp.
    j=nullspace(nullspace(hp,256)+nullspace(anchor,256),256)
    assert len(j)==61 and all(contains(hp,x) and contains(anchor,x) for x in j)
    assert len(echelon(hp+anchor))==93 and all(contains(hq,x) for x in hp+anchor)
    assert contains(j,(1<<256)-1) and all(x.bit_count()%4==0 for x in j)
    cj,ca=cyclic_description(j),cyclic_description(anchor)
    factor,rem=binary_poly_divmod(int(cj['generator_hex'],16),int(ca['generator_hex'],16))
    assert rem==0 and factor==0x169
    seen,x=set(),1
    while x not in seen:
        seen.add(x)
        x=multiply_by_x_mod(x,factor)
    assert x==1 and seen==set(range(1,256))
    model,scales=previous_build()
    variables=model['metadata']['variables']
    variables += [f't_{w}' for w in range(0,129,4)]
    for w in range(0,129,4):
        scales[f't_{w}']=max(1,math.comb(256,w)//(1<<194))
    a=model['metadata']['Lhull_anchor']
    known={int(w):int(v) for w,v in a['half_spectrum'].items()}
    known_dual={int(w):int(v) for w,v in a['dual_half_spectrum'].items()}
    kt=kraw_table()

    def add(name,terms,sense='ge',rhs=0):
        coefficients={}
        for scalar,values in terms:
            for n,v in values.items():
                coefficients[n]=coefficients.get(n,0)+scalar*v
        model['constraints'].append(dict(name=name,coeffs={n:str(v) for n,v in coefficients.items() if v},
                                        sense=sense,rhs=str(rhs)))

    def transform(prefix,degree):
        return {f'{prefix}_{w}':kt[degree][w]+(kt[degree][256-w] if w!=128 else 0)
                for w in range(0,129,4)}

    add('J_mass',[(1,transform('t',0))],'eq',1<<61)
    for w in range(0,129,4):
        t={f't_{w}':1}
        add('J_le_HP_'+str(w),[(1,{f'r_{w}':1}),(-1,t)])
        add('J_le_anchor_'+str(w),[(1,t)],'le',known[w])
        # The common nonzero anchor/J coset is contained in the matching HQ/HP coset.
        add('anchor_J_coset_le_s_'+str(w),[(1,t),(255,{f's_{w}':1})],'ge',known[w])
        if known[w]==0 or w==0:
            add('J_support_'+str(w),[(1,t)],'eq',int(w==0))
    for degree in range(0,129,2):
        t,r,s=(transform(a,degree) for a in ('t','r','s'))
        # Since HP+anchor=HQ, HPdual intersect anchor-dual=HQdual. Their
        # union is contained in Jdual, giving this inclusion-exclusion bound.
        add('Jdual_union_lower_'+str(degree),
            [(1<<32,t),(-255,r),(255,s)],'ge',(1<<93)*known_dual[degree])
        # Root-based BCH bound for Jdual can be added separately after audit.
        # J is doubly even, hence self-orthogonal and J subset Jdual.
        if degree%4==0:
            add('J_self_orthogonal_'+str(degree),[(1,t),(-(1<<61),{f't_{degree}':1})])
    model['metadata']['J_intersection']={
        'dimension':61,'basis_hex':[hex(x) for x in j],
        'cyclic_description':cj,'anchor_cyclic_description':ca,
        'anchor_quotient_factor_hex':hex(factor),'nonzero_shift_orbit':255,
        'HP_plus_anchor_equals_HQ':True,
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/audit_bch_hulls.py',ROOT/'code/prepare_bch_hull_anchor_probe.py')}}
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
