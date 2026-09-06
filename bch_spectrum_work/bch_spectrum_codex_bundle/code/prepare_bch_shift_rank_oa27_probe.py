"""Strength-27 Q/coset moments from the independent distance-28 certificate."""
import json
from pathlib import Path
from prepare_bch_shift_rank_probe import build as previous_build
from certify_bch_shift_rank_threshold import build as rank_build
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_probe import export
from bch_hull_cut_iteration import warm_basis,solve
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/shift_rank_oa27_hull_probe'
PREVIOUS=ROOT/'generated/shift_rank_hull_probe'


def build():
    proof=rank_build(39,'shift_rank_q30_complete')
    source=ROOT/'generated/bch256_shift_rank_q28.json'
    assert proof==json.loads(source.read_text()) and proof['extended_dual_distance_lower']>=28
    model,scales=previous_build()
    kt=kraw_table()
    for j in (24,26):
        for a in ('q','h'):
            model['constraints'].append(dict(name=f'CERTIFIED_OA27_{a}_{j}',
                coeffs={f'{a}_{w}':str(kt[j][w]+(kt[j][256-w] if w!=128 else 0)) for w in range(0,129,2)},
                sense='eq',rhs='0'))
    model['metadata']['Q_OA27_independent_proof']={
        'dual_distance_lower':28,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_shift_rank_threshold.py',source)}}
    return model,scales


if __name__=='__main__':
    model,scales=build()
    lp,count=export(model,scales)
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'scales.json',scales)
    old=json.loads((PREVIOUS/'model.json').read_text())
    basis,mapping=warm_basis(old,scales,model,scales,PREVIOUS/'h_38.bas')
    for name,text in (('h_38.lp',lp),('warm.bas',basis)):
        with (FOLDER/name).open('x') as f:
            f.write(text)
    write_new(FOLDER/'warm_mapping.json',mapping)
    print(json.dumps(dict(variables=len(scales),rows=count)),flush=True)
    print(solve(FOLDER,90,True),flush=True)
