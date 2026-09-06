"""Exact rational shell objectives on the fully audited BCH sandwich model.

The floating LP screen is ill-conditioned; only rational primal/dual
checks in this file produce new shell caps. Outputs are local and write-once.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq
import bridge as base
import occupation_two
from prepare_bch_shift_rank_oa29_probe import build as model_build
from prepare_bch_hull_probe import export
from bch_hull_cut_iteration import solve
from verify_scaled_rational_solution import normalized_rows,parse_assignments


def fq(value):
    return fmpq(value.numerator,value.denominator)


def audit(folder,w,model,scales):
    text=(folder/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    variables=model['metadata']['variables'];rows=normalized_rows(model,scales)
    primal=parse_assignments(text,'VARS:','REDUCED COST:');dual=parse_assignments(text,'PI:','SLACK:')
    assert set(primal)<=set(variables) and set(dual)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={v:fq(primal.get(v,F(0))) for v in variables};assert all(v>=0 for v in x.values())
    combined={v:fmpq(0) for v in variables};bound=fmpq(0)
    for i,row in enumerate(rows,1):
        lhs=sum((fq(a)*x[v] for v,a in row['coeffs'].items()),fmpq(0));rhs=fq(row['rhs'])
        assert {'eq':lhs==rhs,'le':lhs<=rhs,'ge':lhs>=rhs}[row['sense']],row['name']
        price=fq(dual.get(f'c{i}',F(0)))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        bound+=price*rhs
        for v,a in row['coeffs'].items():combined[v]+=price*fq(a)
    target={f'q_{w}':1,f'h_{w}':31}
    assert all(combined[v]>=target.get(v,0) for v in variables)
    assert bound==sum((x[v]*a for v,a in target.items()),fmpq(0))
    assert scales[f'q_{w}']==scales[f'h_{w}']
    physical=F(int(bound.numerator),int(bound.denominator))*scales[f'q_{w}']
    return dict(weight=w,cap=math.floor(physical),physical_upper=base.encode(physical),
                rows_checked=len(rows),variables_checked=len(variables),rational_primal_dual_checks_passed=True)


def run(weights,verify=False):
    reference=base.BCH/'generated/shift_rank_oa29_joint'
    model,scales=model_build()
    assert model==base.read(reference/'model.json') and scales==base.read(reference/'scales.json')
    original,count=export(model,scales);assert original.count(' obj: h_38\n')==1
    for w in weights:
        folder=base.HERE/'generated'/f'joint_shell_exact_w{w}'
        lp=original.replace(' obj: h_38\n',f' obj: q_{w} + 31 h_{w}\n')
        if not verify:
            folder.mkdir(exist_ok=False)
            with (folder/'h_38.lp').open('x') as stream:stream.write(lp)
            with (folder/'warm.bas').open('x') as stream:stream.write((reference/'h_38.bas').read_text())
            print('Solving exact shell',w,flush=True)
            solve(folder,120,True)
        assert (folder/'h_38.lp').read_text()==lp
        result=audit(folder,w,model,scales)
        paths=[Path(__file__),folder/'h_38.lp',folder/'h_38.sol',folder/'warm.bas']
        result['local_sha256']={str(p.relative_to(base.HERE)):base.sha(p) for p in paths}
        result['outer_sha256']={str(p.relative_to(base.BCH)):base.sha(p) for p in
            (reference/'model.json',reference/'scales.json',reference/'audit.json',
             base.BCH/'code/prepare_bch_shift_rank_oa29_probe.py',base.BCH/'code/verify_scaled_rational_solution.py')}
        if verify:assert result==base.read(folder/'cap.json')
        else:base.write_new(folder/'cap.json',result)
        print('Exact shell',w,'cap bits',math.log2(result['cap']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--weights',nargs='+',type=int,required=True);p.add_argument('--verify',action='store_true')
    a=p.parse_args();run(a.weights,a.verify)
