"""Add certified local moments, eliminate fixed variables, select equality rows.

Modular row selection is used only to remove constraints. No rational
equivalence is claimed; any removed row can be checked on a final primal.
"""
import copy
import json
import math
import sys
import textwrap
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
from flint import nmod_mat
from prepare_bch_split38_probe import kraw,export
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_local_probe'


def build():
    model=json.loads((ROOT/'generated/split38_probe/model.json').read_text())
    evidence=json.loads((ROOT/'generated/bch256_wambach_shortening.json').read_text())
    assert evidence['projection']['projection_rank']==38
    assert evidence['checks']['shortened_dual_minimum_at_least']==7
    rows=[r for r in model['constraints'] if not r['name'].startswith('projection_fiber_')]
    kt=kraw(218)
    for i in range(20):
        for degree in range(7):
            coefficients={name:sum(kt[b][degree] for a,b in orbit if a==i)
                          for name,orbit in model['split_orbits'].items()}
            coefficients={n:str(v) for n,v in coefficients.items() if v}
            rhs=math.comb(38,i)*(1<<85) if degree==0 else 0
            if coefficients:
                rows.append(dict(name=f'local_moment_{i}_{degree}',coeffs=coefficients,sense='eq',rhs=str(rhs)))
            else:
                assert rhs==0
    fixed={}
    while True:
        updated=[]
        changed=False
        for row in rows:
            coefficients={n:int(v) for n,v in row['coeffs'].items() if n not in fixed}
            rhs=Fraction(int(row['rhs']))-sum((int(v)*fixed[n] for n,v in row['coeffs'].items() if n in fixed),Fraction(0))
            assert rhs.denominator==1
            if not coefficients:
                assert {'eq':rhs==0,'ge':rhs<=0,'le':rhs>=0}[row['sense']],row['name']
                continue
            if row['sense']=='eq' and len(coefficients)==1:
                name,coeff=next(iter(coefficients.items()))
                value=rhs/coeff
                assert value>=0 and value.denominator==1
                fixed[name]=value
                changed=True
                continue
            updated.append(dict(name=row['name'],coeffs={n:str(v) for n,v in coefficients.items()},sense=row['sense'],rhs=str(rhs.numerator)))
        rows=updated
        if not changed:
            break
    variables=[n for n in model['variables'] if n not in fixed]
    equations=[r for r in rows if r['sense']=='eq']
    index={n:i for i,n in enumerate(variables)}
    prime=2147483647
    transposed=[[0]*len(equations) for _ in range(len(variables)+1)]
    for j,row in enumerate(equations):
        for name,value in row['coeffs'].items():
            transposed[index[name]][j]=int(value)%prime
        transposed[-1][j]=int(row['rhs'])%prime
    reduced,rank=nmod_mat(transposed,prime).rref()
    selected=[]
    for i in range(rank):
        selected.append(next(j for j in range(len(equations)) if reduced[i,j]))
    keep={equations[i]['name'] for i in selected}
    original_rows=copy.deepcopy(rows)
    rows=[r for r in rows if r['sense']!='eq' or r['name'] in keep]
    return dict(classification='Code-specific split LP with necessary local moments; no cap certified yet',
                variables=variables,scales={n:model['scales'][n] for n in variables},
                split_orbits=model['split_orbits'],constraints=rows,objective=model['objective'],
                fixed_physical_variables={n:str(v) for n,v in fixed.items()},
                equality_selection=dict(prime=prime,original=len(equations),retained=rank,
                                        interpretation='Deleting rows only weakens the relaxation; no exact equivalence claim'),
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                               (Path(__file__),ROOT/'generated/split38_probe/model.json',ROOT/'generated/bch256_wambach_shortening.json')}),original_rows


if __name__=='__main__':
    FOLDER.mkdir(exist_ok=False)
    model,all_rows=build()
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'all_rows_after_substitution.json',all_rows)
    export(model,FOLDER/'h38_unwrapped.lp')
    raw=(FOLDER/'h38_unwrapped.lp').read_text()
    wrapped='\n'.join('\n'.join(textwrap.wrap(line,width=2000,break_long_words=False,break_on_hyphens=False,
                                             subsequent_indent=' ')) for line in raw.splitlines())+'\n'
    assert raw.split()==wrapped.split()
    with (FOLDER/'h38.lp').open('x',encoding='ascii') as stream:
        stream.write(wrapped)
    print(json.dumps(dict(variables=len(model['variables']),rows=len(model['constraints']),
                          fixed_variables=len(model['fixed_physical_variables']),
                          equality_selection=model['equality_selection']),indent=2))
