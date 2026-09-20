"""Fix the old exact primal and solve only for its J-spectrum extension."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_intersection_probe import build as full_build
from verify_scaled_rational_solution import parse_assignments
from export_scaled_rational_lp import expression
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_fixed_extension'
PREVIOUS=ROOT/'generated/oa21_hull_reciprocal_probe'


def build():
    model,scales=full_build()
    oldmodel=json.loads((PREVIOUS/'model.json').read_text())
    oldscales=json.loads((PREVIOUS/'scales.json').read_text())
    text=(PREVIOUS/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    raw=parse_assignments(text,'VARS:','REDUCED COST:')
    fixed={n:raw.get(n,Fraction(0))*oldscales[n] for n in oldmodel['metadata']['variables']}
    for row in model['constraints']:
        if row['name'].startswith('J_support_'):
            n,c=next(iter(row['coeffs'].items()))
            assert int(c)==1
            fixed[n]=Fraction(int(row['rhs']))
    variables=[n for n in model['metadata']['variables'] if n not in fixed]
    constraints=[]
    for row in model['constraints']:
        rhs=Fraction(int(row['rhs']))-sum((int(c)*fixed[n] for n,c in row['coeffs'].items() if n in fixed),Fraction(0))
        coefficients={n:int(c)*scales[n] for n,c in row['coeffs'].items() if n not in fixed and int(c)}
        if not coefficients:
            assert {'eq':rhs==0,'ge':rhs<=0,'le':rhs>=0}[row['sense']],row['name']
            continue
        norm=max(Fraction(1),abs(rhs),*(abs(c) for c in coefficients.values()))
        constraints.append(dict(name=row['name'],coeffs={n:str(Fraction(c)/norm) for n,c in coefficients.items()},
                                sense=row['sense'],rhs=str(rhs/norm)))
    return dict(variables=variables,scales={n:scales[n] for n in variables},
                fixed_physical={n:str(v) for n,v in fixed.items()},constraints=constraints,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                    (Path(__file__),ROOT/'code/prepare_bch_hull_intersection_probe.py',
                     PREVIOUS/'model.json',PREVIOUS/'scales.json',PREVIOUS/'h_38.sol',PREVIOUS/'audit.json')})


if __name__=='__main__':
    model=build()
    lines=['Maximize',' obj: t_64','Subject To']
    for i,row in enumerate(model['constraints'],1):
        pieces=[expression([(Fraction(v),n)]) for n,v in row['coeffs'].items()]
        chunks=[' '.join(pieces[j:j+4]) for j in range(0,len(pieces),4)]
        chunks[0]=f' c{i}: '+chunks[0]
        chunks[-1]+=' '+{'eq':'=','ge':'>=','le':'<='}[row['sense']]+' '+row['rhs']
        lines+=chunks
    lines+=['Bounds']+[f' 0 <= {n}' for n in model['variables']]+['End']
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    with (FOLDER/'h_38.lp').open('x') as f:
        f.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(variables=len(model['variables']),rows=len(model['constraints']))))
