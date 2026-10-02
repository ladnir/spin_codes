"""Retune the actual small-step fixed-input bound at one dense bottleneck."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import math
from scipy.optimize import minimize_scalar
from flint import arb,ctx
import model
import target49_dense
import retune_fixed_input


def run(a):
    if a.output.exists(): raise FileExistsError(a.output)
    ctx.prec=256
    c=target49_dense.Checker(a.t,a.s)
    q,v=a.q,F(a.coordinate)
    best=None
    for index in (-120,-100,-80,-60,-40,-24,-12,0,12):
        def evaluate(u):
            nonlocal best
            witness=retune_fixed_input.proposal(c,q,v,index,float(u))
            value=c.fixed_input_bound(q,q,v,v,witness)
            score=float(value)
            if best is None or score<best[0]: best=(score,witness)
            return score
        grid=[-2,-1,-.5,0,.5,1,2]
        scores=[evaluate(u) for u in grid]
        i=min(range(len(grid)),key=lambda i:scores[i])
        minimize_scalar(evaluate,bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),
                        method='bounded',options=dict(maxiter=18,xatol=1e-4))
        print('retune',index,'best margin',-best[0]/math.log(2),flush=True)
    upper=c.bound(q,q,v,v,best[1]).exp()
    result=dict(status='K16_TARGET49_RETUNED_DENSE_POINT',instance=c.engine.identity(),
        q=q,coordinate=str(v),witness=best[1],upper=model.base.pack(upper),
        margin_bits=float(-upper.log()/arb(2).log()),full_distance_proved=False,
        source_sha256=model.base.sources())
    for module in (Path(__file__),Path(target49_dense.__file__),Path(retune_fixed_input.__file__)):
        result['source_sha256'][module.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(module)
    model.base.base.write_new(a.output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--t',type=int,default=32)
    p.add_argument('--s',type=int,default=15)
    p.add_argument('--q',type=int,default=218)
    p.add_argument('--coordinate',default='21/128')
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
