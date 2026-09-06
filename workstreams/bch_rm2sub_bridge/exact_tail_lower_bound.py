"""Minimize the actual BCH low-weight tail in the retained rational model.

Maximize its negative. A dual upper bound on the negative objective is a
lower bound on the tail; it is never fed to the shell-upper-cap loader.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq
import bridge as base
import exact_joint_shell_caps as exact
from export_scaled_rational_lp import expression


def run(verify=False):
    reference=base.BCH/'generated/shift_rank_oa29_joint'
    model,scales=exact.model_build()
    assert model==base.read(reference/'model.json') and scales==base.read(reference/'scales.json')
    original,_=exact.export(model,scales);folder=base.HERE/'generated/joint_tail_lower_80'
    common=scales['q_80']
    target={f'{a}_{w}':-F(mult*scales[f'{a}_{w}'],common) for w in range(38,81,2) for a,mult in [('q',1),('h',31)]}
    lp=original.replace(' obj: h_38\n',' obj: '+expression([(value,name) for name,value in target.items()])+'\n')
    if not verify:
        folder.mkdir(exist_ok=False)
        with (folder/'h_38.lp').open('x') as stream:stream.write(lp)
        warm=base.HERE/'generated/joint_shell_chain01_w80/h_38.bas'
        with (folder/'warm.bas').open('x') as stream:stream.write(warm.read_text())
        print('Solving exact lower bound on BCH weights 38..80',flush=True)
        exact.solve(folder,300,True)
    assert (folder/'h_38.lp').read_text()==lp
    solution=(folder/'h_38.sol').read_text();assert 'status = OPTIMAL' in solution
    variables=model['metadata']['variables'];rows=exact.normalized_rows(model,scales)
    primal=exact.parse_assignments(solution,'VARS:','REDUCED COST:');dual=exact.parse_assignments(solution,'PI:','SLACK:')
    assert set(primal)<=set(variables) and set(dual)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={v:exact.fq(primal.get(v,F(0))) for v in variables};assert all(v>=0 for v in x.values())
    combined={v:fmpq(0) for v in variables};bound=fmpq(0)
    for i,row in enumerate(rows,1):
        lhs=sum((exact.fq(a)*x[v] for v,a in row['coeffs'].items()),fmpq(0));rhs=exact.fq(row['rhs'])
        assert {'eq':lhs==rhs,'le':lhs<=rhs,'ge':lhs>=rhs}[row['sense']],row['name']
        price=exact.fq(dual.get(f'c{i}',F(0)))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        bound+=price*rhs
        for v,a in row['coeffs'].items():combined[v]+=price*exact.fq(a)
    assert all(combined[v]>=exact.fq(target.get(v,F(0))) for v in variables)
    assert bound==sum((x[v]*exact.fq(a) for v,a in target.items()),fmpq(0))
    lower=-F(int(bound.numerator),int(bound.denominator))*common;integer=max(0,math.ceil(lower))
    result=dict(status='EXACT_BCH_TAIL_LOWER',weights=list(range(38,81,2)),lower=integer,
        rational_lower=base.encode(lower),rational_primal_dual_checks_passed=True,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),Path(exact.__file__),
            folder/'h_38.lp',folder/'h_38.sol',folder/'warm.bas']},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in [reference/'model.json',reference/'scales.json',reference/'audit.json',
            base.BCH/'code/prepare_bch_shift_rank_oa29_probe.py',base.BCH/'code/verify_scaled_rational_solution.py']})
    if verify:assert result==base.read(folder/'lower.json')
    else:base.write_new(folder/'lower.json',result)
    print('Exact BCH tail lower bound',integer,'log2',math.log2(integer) if integer else None,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
