"""Improve an existing exact BCH dual witness using small local moment LPs.

Floating optimization proposes degree-six polynomials only. Their pointwise
inequalities and the final global upper bound are checked with exact fractions.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import numpy as np
from scipy.optimize import linprog
from flint import fmpq,fmpq_mat
from audit_bch_q1_full_arb import encode,decode
from prepare_bch_split38_probe import kraw
from verify_scaled_rational_solution import normalized_rows,parse_assignments
from export_scaled_rational_lp import variable_scales
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
OUTPUT=GEN/'bch256_local_slack_certificate.json'


def residuals():
    folder=GEN/'oa21_closure_probe'
    model=json.loads((folder/'model.json').read_text())
    meta=json.loads((folder/'h_38.json').read_text())
    scale=meta['objective_physical_multiplier']
    scales=variable_scales(model)
    prices=parse_assignments((folder/'h_38_dual.sol').read_text(),'PI:','SLACK:')
    combined={name:Fraction(0) for name in model['metadata']['variables']}
    upper=Fraction(0)
    for i,row in enumerate(normalized_rows(model,scales),1):
        price=prices.get(f'c{i}',Fraction(0))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        upper+=price*row['rhs']*scale
        for name,value in row['coeffs'].items():
            combined[name]+=price*value*scale/scales[name]
    combined['h_38']-=1
    assert all(v>=0 for v in combined.values())
    saved=json.loads((folder/'audit.json').read_text())
    assert upper==decode(saved['shells'][0]['physical_optimum'])
    return combined,upper


def zero_primal(basis,costs):
    allowed=[j for j,c in enumerate(costs) if c==0]
    if not allowed:
        return None
    result=linprog(np.zeros(len(allowed)),A_eq=np.array([[float(row[j]) for j in allowed] for row in basis]),
                   b_eq=[1]+[0]*6,bounds=(0,None),method='highs')
    if not result.success:
        return None
    selected=[allowed[i] for i,x in enumerate(result.x) if x>1e-9]
    if not selected:
        return None
    matrix=fmpq_mat([[fmpq(basis[t][j].numerator,basis[t][j].denominator) for j in selected] for t in range(7)])
    transposed,rank=matrix.transpose().rref()
    if rank!=len(selected):
        return None
    pivot_rows=[next(t for t in range(7) if transposed[i,t]) for i in range(rank)]
    square=fmpq_mat([[matrix[t,j] for j in range(rank)] for t in pivot_rows])
    rhs=fmpq_mat([[int(t==0)] for t in pivot_rows])
    solution=square.solve(rhs)
    values=[Fraction(int(solution[i,0].numerator),int(solution[i,0].denominator)) for i in range(rank)]
    if any(v<0 for v in values):
        return None
    if any(sum((basis[t][j]*x for j,x in zip(selected,values)),Fraction(0))!=int(t==0) for t in range(7)):
        return None
    assert all(costs[j]==0 for j in selected)
    return [dict(node_index=j,mass=encode(x)) for j,x in zip(selected,values)]


def build():
    residual,old=residuals()
    kt=kraw(218)
    blocks=[]
    total=Fraction(0)
    def symmetric_cost(prefix,w):
        return residual[f'{prefix}_{min(w,256-w)}']/(1 if w==128 else 2)
    for i in range(20):
        nodes=[j for j in range(219) if (i+j)%2==0
               and (i+j in (0,256) or 40<=i+j<=216) and 38<=38-i+j<=218]
        costs=[symmetric_cost('q',i+j)+symmetric_cost('h',38-i+j) for j in nodes]
        basis=[[Fraction(kt[j][t],math.comb(218,t)) for j in nodes] for t in range(7)]
        best=[Fraction(0)]*7
        lower=Fraction(0)
        proposals=0
        for exponent in (-100,-90,-80,-70,-60,-50,-40):
            unit=Fraction(1,1<<(-exponent))
            objective=np.array([float(min(Fraction(1),c/unit)) for c in costs])
            result=linprog(objective,A_eq=np.array([[float(v) for v in row] for row in basis]),
                           b_eq=[1]+[0]*6,bounds=(0,None),method='highs')
            if not result.success:
                continue
            proposals+=1
            coefficients=[Fraction.from_float(float(v))*unit for v in result.eqlin.marginals]
            violation=max([Fraction(0)]+[sum((coefficients[t]*basis[t][j] for t in range(7)),Fraction(0))-costs[j]
                                        for j in range(len(nodes))])
            coefficients[0]-=violation
            if coefficients[0]>lower:
                best,lower=coefficients,coefficients[0]
        assert all(sum((best[t]*basis[t][j] for t in range(7)),Fraction(0))<=costs[j] for j in range(len(nodes)))
        multiplicity=1 if i==19 else 2
        mass=math.comb(38,i)*(1<<85)
        contribution=multiplicity*mass*lower
        total+=contribution
        primal=zero_primal(basis,costs)
        blocks.append(dict(inside_weight=i,outside_nodes=nodes,complement_multiplicity=multiplicity,
                           fiber_mass=str(mass),normalized_kraw_polynomial_coefficients=[encode(v) for v in best],
                           certified_slack_lower=encode(contribution),floating_proposals_checked=proposals,
                           exact_zero_cost_primal=primal))
        print('LOCAL',i,'slack',float(contribution),'zero_cost_witness',primal is not None,flush=True)
    assert 0<=total<=old
    upper=old-total
    cap=31*(upper.numerator//upper.denominator)
    paths=[Path(__file__),ROOT/'code/verify_scaled_rational_solution.py',ROOT/'code/export_scaled_rational_lp.py',
           ROOT/'code/export_lp.py',ROOT/'code/prepare_bch_split38_probe.py',GEN/'oa21_closure_probe/model.json',
           GEN/'oa21_closure_probe/h_38.json',GEN/'oa21_closure_probe/h_38_dual.sol',GEN/'oa21_closure_probe/audit.json',
           GEN/'bch256_wambach_shortening.json']
    return dict(classification='Exact local-polynomial correction of the existing OA21 dual witness',
                old_h38_upper=encode(old),certified_total_slack=encode(total),new_h38_upper=encode(upper),
                A38_cap_without_lattice_rounding=cap,A38_cap_at_most_10_to_13=cap<=10**13,
                all_local_polynomials_checked_exactly=True,
                blocks_with_exact_zero_cost_primal=sum(b['exact_zero_cost_primal'] is not None for b in blocks),
                blocks=blocks,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:result[k] for k in ('A38_cap_without_lattice_rounding','A38_cap_at_most_10_to_13',
                                          'blocks_with_exact_zero_cost_primal')},indent=2))
