"""Use the complete exact count A16(HQdual)=16592."""
import json
from pathlib import Path
from certify_bch_hull_dual_weight16 import build as shell_build
from prepare_bch_hull_dual16_probe import build as previous_build
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_probe import export
from bch_hull_cut_iteration import warm_basis,solve
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_dual16_count_probe'
PREVIOUS=ROOT/'generated/oa21_hull_dual16_probe'


def build():
    shell=shell_build()
    assert shell==json.loads((ROOT/'generated/bch256_hull_dual_weight16.json').read_text())
    assert shell['full_Adual16_shell_enumerated'] and shell['HQdual_A16_exact']==16592
    model,scales=previous_build()
    kt=kraw_table()
    coeffs={f'{a}_{w}':str(mult*(kt[16][w]+(kt[16][256-w] if w!=128 else 0)))
            for a,mult in (('r',1),('s',255)) for w in range(0,129,4)}
    model['constraints'].append(dict(name='EXACT_HQdual_weight16_count',coeffs=coeffs,sense='eq',rhs=str((1<<93)*16592)))
    model['metadata']['HQdual16_enumeration']={'exact_shell_count':16592,
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_hull_dual_weight16.py',
             ROOT/'generated/bch256_hull_dual_weight16.json',ROOT/'code/prepare_bch_hull_dual16_probe.py')}}
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
