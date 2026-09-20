"""Add the complete exact zero weight-18 shell of HQdual."""
import json
from pathlib import Path
from certify_bch_hull_affine_shell import build as shell_build
from prepare_bch_hull_dual16_count_probe import build as previous_build
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_probe import export
from bch_hull_cut_iteration import warm_basis,solve
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_dual18_probe'
PREVIOUS=ROOT/'generated/oa21_hull_dual16_count_probe'


def build():
    shell=shell_build(18,'hull_weight18_affine_triples')
    source=ROOT/'generated/bch256_hull_dual_weight18_complete.json'
    assert shell==json.loads(source.read_text())
    assert shell['full_Adual_shell_enumerated'] and shell['HQdual_shell_count_exact']==0
    model,scales=previous_build()
    kt=kraw_table()
    coefficients={f'{a}_{w}':str(mult*(kt[18][w]+(kt[18][256-w] if w!=128 else 0)))
                  for a,mult in (('r',1),('s',255)) for w in range(0,129,4)}
    model['constraints'].append(dict(name='EXACT_HQdual_weight18_zero',coeffs=coefficients,sense='eq',rhs='0'))
    model['metadata']['HQdual18_enumeration']={'exact_shell_count':0,
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_hull_affine_shell.py',source,
             ROOT/'code/prepare_bch_hull_dual16_count_probe.py')}}
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
