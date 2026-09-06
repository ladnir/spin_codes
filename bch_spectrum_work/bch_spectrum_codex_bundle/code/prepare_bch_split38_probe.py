"""Exploratory split enumerator around a weight-38 word in P minus Q.

Only necessary continuous constraints are encoded. No new BCH cap is claimed
without an exact witness and an applicability audit. Saved models are new.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
from export_lp import independent_oa_constraints
from export_scaled_rational_lp import expression,fraction_text
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_probe'


def kraw(n):
    result=[]
    for x in range(n+1):
        col=[1,n-2*x]
        for d in range(1,n):
            value,remainder=divmod((n-2*x)*col[-1]-(n-d+1)*col[-2],d+1)
            assert not remainder
            col.append(value)
        result.append(col[:n+1])
    return result


def build():
    old=json.loads((ROOT/'generated/oa21_closure_probe/model.json').read_text())
    variables=old['metadata']['variables'].copy()
    rows=[row.copy() for row in independent_oa_constraints(old)
          if not row['name'].startswith(('q_nonneg_','h_nonneg_'))]
    splits={}
    for i in range(39):
        for j in range(219):
            if (i+j)%2:
                continue
            if i+j not in (0,256) and not 40<=i+j<=216:
                continue
            if not 38<=38-i+j<=218:
                continue
            pair=(38-i,218-j)
            if (i,j)>pair:
                continue
            splits[f'b_{i}_{j}']=[(i,j),pair]
    variables.extend(splits)
    def add(name,co,sense='eq',rhs=0):
        rows.append(dict(name=name,coeffs={n:str(v) for n,v in co.items() if v},sense=sense,rhs=str(rhs)))
    for prefix in ('q','h'):
        for w in range(0,129,2):
            co={f'{prefix}_{w}':-1}
            for name,orbit in splits.items():
                co[name]=sum(int((i+j if prefix=='q' else 38-i+j)==w) for i,j in orbit)
            add(f'split_{prefix}_{w}',co)
    ka,kb=kraw(38),kraw(218)
    for total in range(2,25,2):
        for u in range(min(38,total)+1):
            v=total-u
            co={name:sum(ka[i][u]*kb[j][v] for i,j in orbit) for name,orbit in splits.items()}
            add(f'split_dual_{u}_{v}',co,'eq' if total<=20 else 'ge')
    # The dual shortened to 38 positions has dimension <=2: eight words of
    # distance >=22 contradict 28*22 > 38*16. Thus projection rank >=36.
    for i in range(20):
        co={name:sum(int(a==i) for a,b in orbit) for name,orbit in splits.items()}
        add(f'projection_fiber_{i}',co,'le',math.comb(38,i)*(1<<87))
    scales={}
    for name in variables:
        if name in splits:
            i,j=splits[name][0]
            scales[name]=max(1,math.comb(38,i)*math.comb(218,j)//(1<<132))
        else:
            scales[name]=max(1,math.comb(256,int(name.split('_')[1]))//(1<<132))
    return dict(classification='Provisional necessary split-enumerator LP, not a certified bound',
                variables=variables,scales=scales,split_orbits=splits,constraints=rows,
                objective=dict(variable='h_38',physical_multiplier=scales['h_38']),
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                               (Path(__file__),ROOT/'generated/oa21_closure_probe/model.json')})


def export(model,path):
    variables,scales=model['variables'],model['scales']
    lines=['Maximize',' obj: h_38','Subject To']
    for row in model['constraints']:
        terms=[(int(v)*scales[n],n) for n,v in row['coeffs'].items()]
        rhs=int(row['rhs'])
        norm=max(1,abs(rhs),*(abs(v) for v,n in terms))
        terms=[(Fraction(v,norm),n) for v,n in terms]
        op={'eq':'=','ge':'>=','le':'<='}[row['sense']]
        lines.append(f' {expression(terms)} {op} {fraction_text(Fraction(rhs,norm))}')
    lines+=['Bounds']+[f' 0 <= {n}' for n in variables]+['End']
    with path.open('x',encoding='ascii') as stream:
        stream.write('\n'.join(lines)+'\n')


if __name__=='__main__':
    FOLDER.mkdir(exist_ok=False)
    model=build()
    write_new(FOLDER/'model.json',model)
    export(model,FOLDER/'h38.lp')
    print(json.dumps(dict(variables=len(model['variables']),rows=len(model['constraints']),
                          nonzeros=sum(len(r['coeffs']) for r in model['constraints'])),indent=2))
