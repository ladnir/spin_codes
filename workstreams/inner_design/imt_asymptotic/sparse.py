"""First-order IMT Collatz witness and finite-alpha diagnostics; not a theorem."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
import screen


def witness(model):
    e = model.engine
    m = e.m
    v = F(1,1024)
    beta_scale,lambda_scale = F(4,5),F(8,5)
    # d_i counts single-position syndromes whose expansion is in shell i.
    d = [sum(e.low[1]['by_weight'].get(w,{}).values()) for w in e.levels]
    assert sum(d)==128
    reset = [beta_scale/2*(F(k,n)+F(128,m))/v for k,n in zip(d,map(int,model.counts))]
    slopes = [r-lambda_scale*w for r,w in zip(reset,e.levels)]
    mean = sum(F(int(n),m)*s for n,s in zip(model.counts,slopes))
    h = [2*(s-mean) for s in slopes]
    assert sum(F(int(n),m)*c for n,c in zip(model.counts,h))==0
    return v,[F(900),*h],mean


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--gamma',type=int,default=96)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    model = screen.Model()
    v,h,mean = witness(model)
    rows = []
    for alpha in np.geomspace(1e-8,1e-4,81):
        beta = .8*alpha
        lam = -math.log1p(-1.6*alpha)
        matrix = model.occupation(beta,lam)
        vector = np.r_[1.,float(v)*(1+alpha*np.array(list(map(float,h))))]
        residual = (matrix @ vector-(1-args.gamma*alpha)*vector)/(alpha*vector)
        exponent = alpha*math.log(2)/2+screen.kl(alpha,1.6*alpha)+model.rate(beta,lam,'occupation')+.11*lam
        rows.append(dict(alpha=float(alpha),max_normalized_collatz_residual=float(max(residual)),
                         random_exponent_per_alpha=exponent/alpha,
                         structured_with_support_per_alpha=exponent/alpha+.01281+math.log(2)/9.75))
    # The reduced impulse envelope deliberately does not subtract the reset
    # mass from survival. This matches the existing nonnegative upper transfer.
    reduced = [[0.,1.],[1/model.m,1.]]
    first = [-128*F(4,5)*(1-v)+args.gamma]
    reset_d = F(4,5)*(F(1,2)+F(64,model.m))/v
    first.append(-h[0]/2-F(8,5)*48+reset_d+args.gamma)
    first.extend([mean+args.gamma]*len(model.levels))
    assert all(x<0 for x in first)
    result = dict(status='IMT_SPARSE_FIRST_ORDER_AND_POINT_SCREEN_NOT_CERTIFICATE',
        gamma=args.gamma,p_scale='8/5',beta_scale='4/5',z='1-(8/5)*alpha',v=str(v),
        affine_live_coordinates=[str(x) for x in h],
        exact_first_order_relative_residuals=[str(x) for x in first],
        limiting_live_generator_diagonal=str(mean),reduced_impulse_upper=reduced,
        mean_live_output_per_bit=str(F(1<<18,model.m)),
        sampled_alpha_range=[1e-8,1e-4],rows=rows,
        all_sampled_collatz_negative=all(r['max_normalized_collatz_residual']<0 for r in rows),
        worst=max(rows,key=lambda r:r['max_normalized_collatz_residual']),
        source_sha256={p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in [Path(__file__),Path(screen.__file__)]},
        limitations=['Exact rational first-order coefficients do not control the remainder on a full interval.',
                     'Floating-point samples do not certify alpha approaching zero.',
                     'The continuum reduction and uniform growing-sparse route error remain proof obligations.'])
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_sha256')}))


if __name__=='__main__':
    main()
