"""Floating LP shell objectives on the retained exact BCH constraint model.

Discovery only. No solver decimal objective is a certified shell cap.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import bridge as base
import occupation_two
import christoffel_caps as caps
from verify_scaled_rational_solution import normalized_rows


def run(weights,tag):
    folder=base.BCH/'generated/shift_rank_oa29_joint'
    model=base.read(folder/'model.json');scales=base.read(folder/'scales.json')
    variables=model['metadata']['variables'];indices={v:i for i,v in enumerate(variables)}
    rows=normalized_rows(model,scales);equal=[];eq_rhs=[];inequal=[];in_rhs=[]
    for row in rows:
        a=np.zeros(len(variables))
        for name,value in row['coeffs'].items():a[indices[name]]=float(value)
        rhs=float(row['rhs'])
        if row['sense']=='eq':equal.append(a);eq_rhs.append(rhs)
        else:
            sign=1 if row['sense']=='le' else -1
            inequal.append(sign*a);in_rhs.append(sign*rhs)
    results=[];old=caps.deterministic_caps()
    for w in weights:
        c=np.zeros(len(variables));c[indices[f'q_{w}']]=-1;c[indices[f'h_{w}']]=-31
        answer=linprog(c,A_ub=inequal,b_ub=in_rhs,A_eq=equal,b_eq=eq_rhs,bounds=(0,None),method='highs',
                       options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'time_limit':60})
        row=dict(weight=w,status=int(answer.status),message=answer.message)
        if answer.success:
            physical=-float(answer.fun)*scales[f'q_{w}']
            row.update(log2_cap_diagnostic=math.log2(physical),old_cap_log2=math.log2(old[w]),
                       equality_prices=answer.eqlin.marginals.tolist(),inequality_prices=answer.ineqlin.marginals.tolist())
        results.append(row)
        print({k:v for k,v in row.items() if 'prices' not in k},flush=True)
    base.write_new(base.HERE/'generated'/f'joint_shell_{tag}_screen.json',dict(status='FLOAT_LP_DIAGNOSTIC_ONLY',rows=results,
        source_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in (Path(__file__),)},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in (folder/'model.json',folder/'scales.json',base.BCH/'code/verify_scaled_rational_solution.py')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--weights',nargs='+',type=int,required=True);p.add_argument('--tag',required=True)
    a=p.parse_args();run(a.weights,a.tag)
