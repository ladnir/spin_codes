"""Exact substitution of fixed J coordinates before a higher-precision solve."""
import json
from pathlib import Path
from prepare_bch_hull_intersection_probe import build as original_build
from prepare_bch_hull_probe import export
from run_higher_endpoint_preflight import write_new
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_intersection_reduced'


def build():
    model,scales=original_build()
    fixed={}
    for row in model['constraints']:
        if row['name'].startswith('J_support_'):
            assert row['sense']=='eq' and len(row['coeffs'])==1
            n,c=next(iter(row['coeffs'].items()))
            assert int(c)==1
            fixed[n]=int(row['rhs'])
    newrows=[]
    for row in model['constraints']:
        rhs=int(row['rhs'])-sum(int(c)*fixed[n] for n,c in row['coeffs'].items() if n in fixed)
        coeffs={n:c for n,c in row['coeffs'].items() if n not in fixed and int(c)}
        if not coeffs:
            assert {'eq':rhs==0,'ge':rhs<=0,'le':rhs>=0}[row['sense']]
            continue
        newrows.append(dict(name=row['name'],coeffs=coeffs,sense=row['sense'],rhs=str(rhs)))
    model['constraints']=newrows
    model['metadata']['variables']=[n for n in model['metadata']['variables'] if n not in fixed]
    model['metadata']['exact_fixed_J_substitution']=fixed
    scales={n:s for n,s in scales.items() if n not in fixed}
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
