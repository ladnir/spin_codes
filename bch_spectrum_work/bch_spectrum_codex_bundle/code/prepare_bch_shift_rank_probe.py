"""Certified stronger moments from independently checked Fourier rank witnesses."""
import json
from pathlib import Path
from prepare_bch_hull_dual18_probe import build as previous_build
from certify_bch_shift_rank import build as rank_build
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_probe import export
from verify_scaled_rational_solution import normalized_rows
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/shift_rank_hull_probe'


def build():
    rank=rank_build()
    source=ROOT/'generated/bch256_shift_rank_dual_distances.json'
    assert rank==json.loads(source.read_text())
    assert rank['Q']['extended_dual_distance_lower']>=24 and rank['P']['extended_dual_distance_lower']>=26
    model,scales=previous_build()
    kt=kraw_table()
    def coefficients(degree,multipliers):
        return {f'{a}_{w}':str(mult*(kt[degree][w]+(kt[degree][256-w] if w!=128 else 0)))
                for a,mult in multipliers for w in range(0,129,2)}
    for prefix in ('q','h'):
        model['constraints'].append(dict(name='CERTIFIED_ShiftRank_Qdual24_'+prefix,
            coeffs=coefficients(22,((prefix,1),)),sense='eq',rhs='0'))
    model['constraints'].append(dict(name='CERTIFIED_ShiftRank_Pdual26',
        coeffs=coefficients(24,(('q',1),('h',255))),sense='eq',rhs='0'))
    model['metadata']['independent_dual_rank_proof']={
        'Qdual_distance_lower':24,'Pdual_distance_lower':26,'literature_table_used':False,
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_shift_rank.py',source)}}
    return model,scales


if __name__=='__main__':
    model,scales=build()
    lp,count=export(model,scales)
    prior=ROOT/'generated/schaubplus_provisional'
    old=json.loads((prior/'model.json').read_text())
    oldscales=json.loads((prior/'scales.json').read_text())
    assert scales==oldscales
    oldrows=normalized_rows(old,scales); newrows=normalized_rows(model,scales)
    assert len(oldrows)==len(newrows)
    for a,b in zip(oldrows,newrows):
        assert all(a[k]==b[k] for k in ('coeffs','sense','rhs'))
    assert lp==(prior/'h_38.lp').read_text()
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'scales.json',scales)
    for name,text in [('h_38.lp',lp)]+[(name,(prior/name).read_text()) for name in ('h_38.sol','h_38.bas')]:
        with (FOLDER/name).open('x') as f:
            f.write(text)
    write_new(FOLDER/'witness_reuse.json',dict(exact_exported_LP_identity_verified=True,
        new_optimizer_run=False,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),prior/'h_38.lp',prior/'h_38.sol',prior/'h_38.bas')}))
    print(json.dumps(dict(variables=len(scales),rows=count,exact_previous_witness_reused=True)))
