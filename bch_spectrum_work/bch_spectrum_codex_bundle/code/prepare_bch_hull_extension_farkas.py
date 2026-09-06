"""Find an explicit Farkas witness for the fixed old hull primal."""
import json
from fractions import Fraction
from pathlib import Path
from export_scaled_rational_lp import expression
from prepare_bch_hull_fixed_extension import build as extension_build
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_extension_farkas'


def build():
    source=extension_build()
    assert source==json.loads((ROOT/'generated/oa21_hull_fixed_extension/model.json').read_text())
    variables=[]
    terms=[]
    for i,row in enumerate(source['constraints']):
        for sign in ((1,-1) if row['sense']=='eq' else ((1,) if row['sense']=='ge' else (-1,))):
            n=f'f_{len(variables)}'
            variables.append(n)
            terms.append(dict(variable=n,row_index=i,sign=sign))
    objective={item['variable']:str(item['sign']*Fraction(source['constraints'][item['row_index']]['rhs'])) for item in terms}
    constraints=[]
    for n in source['variables']:
        coeffs={item['variable']:str(item['sign']*Fraction(source['constraints'][item['row_index']]['coeffs'].get(n,'0')))
                for item in terms}
        constraints.append(dict(name=n,coeffs=coeffs,sense='le',rhs='0'))
    constraints.append(dict(name='multiplier_mass',coeffs={n:'1' for n in variables},sense='le',rhs='1'))
    return dict(variables=variables,objective=objective,constraints=constraints,terms=terms,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                    (Path(__file__),ROOT/'generated/oa21_hull_fixed_extension/model.json')})


def wrapped(coeffs):
    pieces=[expression([(Fraction(v),n)]) for n,v in coeffs.items() if Fraction(v)]
    return [' '.join(pieces[j:j+4]) for j in range(0,len(pieces),4)] or ['0 f_0']


if __name__=='__main__':
    model=build()
    lines=['Maximize']+wrapped(model['objective'])+['Subject To']
    for i,row in enumerate(model['constraints'],1):
        chunks=wrapped(row['coeffs'])
        chunks[0]=f' c{i}: '+chunks[0]
        chunks[-1]+=' <= '+row['rhs']
        lines+=chunks
    lines+=['Bounds']+[f' 0 <= {n}' for n in model['variables']]+['End']
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    with (FOLDER/'h_38.lp').open('x') as f:
        f.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(variables=len(model['variables']),rows=len(model['constraints']))))
