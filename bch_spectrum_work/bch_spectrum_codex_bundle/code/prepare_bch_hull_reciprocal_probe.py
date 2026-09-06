"""Add the reciprocal coset containment HQ/HP into Qdual/Pdual."""
import json
from pathlib import Path
from audit_bch_hulls import rows,hull,nullspace,echelon,contains
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_anchor_probe import build as previous_build
from prepare_bch_hull_probe import export
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_reciprocal_probe'


def build():
    p,q=rows(37,131),rows(39,123)
    hp,hq=hull(p),hull(q)
    pd,qd=nullspace(p,256),nullspace(q,256)
    assert all(contains(qd,x) for x in hq)
    assert len(hq)+len(pd)-len(echelon(hq+pd))==85
    assert all(contains(hq,x) and contains(pd,x) for x in hp)
    # HQ intersect Pdual = HP, hence HQ\HP subset Qdual\Pdual.
    model,scales=previous_build()
    kt=kraw_table()
    for j in range(0,129,4):
        coeffs={}
        for prefix,mult in (('q',1),('h',-1)):
            for w in range(0,129,2):
                coeffs[f'{prefix}_{w}']=str(mult*(kt[j][w]+(kt[j][256-w] if w!=128 else 0)))
        coeffs[f's_{j}']=str(-(1<<131))
        model['constraints'].append(dict(name='s_le_Qdual_Pdual_coset_'+str(j),
                                        coeffs=coeffs,sense='ge',rhs='0'))
    model['metadata']['reciprocal_coset']={
        'HQ_intersect_Pdual_dimension':85,'intersection_equals_HP':True,
        'source_sha256':{str(path.relative_to(ROOT)):sha(path) for path in
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
