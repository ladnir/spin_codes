"""Partial paired-map screen for independently sampled transvections.

Each transvection samples u uniformly nonzero and v uniformly in u-perp,
including v=0, then applies R=I+u*v^T. It is invertible since v^T*u=0.
For fixed nonzero a its exact law is half the atom at a plus half uniform
on all nonzero vectors. Independent compositions have lazy weight 2^-r.
This is distinct from an unconstrained rank-one update, which can be singular.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import screen_refresh_stride as stride
from flint import arb,arb_mat,ctx
from packet_outer_geometry_proposal import estimate,upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


def local_family(physical,tilt,updates):
    if updates not in (1,2,3):raise ValueError('one through three transvections required')
    data=dict(physical,updates=updates)
    z=(-stride.screen.aq(Q(tilt))).exp()
    local=stride.lazy.outward_at_z(data,z)
    local=stride.lazy.refine_zero(data,local,z)
    local=stride.screen.q1.kernel_t64.kernel_birth_density.refine_local(data,local,z,Q(1,2))
    # Exact empty-input split, including the tilted uniform source.
    alpha=arb(2)**(-updates);beta=1-alpha
    levels=list(map(int,data['birth_class_levels']))
    n=local[0].nrows();empty=arb_mat(n,n);empty[0,0]=1
    empty[1,1]=alpha*z**min(levels);empty[1,2]=beta*z**min(levels)
    total=arb(0)
    for i,(level,count) in enumerate(zip(levels,data['birth_class_counts']),3):
        mass=stride.screen.up(stride.screen.aq(Q(int(count),65535))*z**level)
        total+=mass
        empty[2,i]=stride.screen.up(alpha*mass)
        empty[i,i]=alpha*z**level;empty[i,2]=beta*z**level
    empty[2,2]=stride.screen.up(beta*total);local[0]=empty
    return stride.screen.q1.kernel_t64.convolve(local)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--updates',type=int,choices=(1,2,3),required=True)
    parser.add_argument('--tilts',nargs='+',required=True)
    parser.add_argument('--max-q',type=int,default=128)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('fresh output required')
    ctx.prec=192
    data,record=stride.maps.prepare()
    physical=stride.screen.q1.kernel_t64.authenticate(data)
    _,counts=stride.screen.tail.uniform_envelope(16,8,4)
    one_best=[arb(1)]*65;best={};choices={};trials=[]
    for tilt in args.tilts:
        local=local_family(physical,tilt,args.updates)
        regional=stride.screen.q1.placement(local,epochs=16,windows=32,rounding=stride.screen.q1.rounded,maximum_groups=1)
        factor=(stride.screen.aq(Q(tilt))*13107).exp()
        moments=stride.screen.q1.support_moments(regional[0],regional[1],64)
        one_best=[min(a,stride.screen.up(factor*b)) for a,b in zip(one_best,moments)]
        one=stride.screen.up(512*sum((stride.screen.aq(c)*m for c,m in zip(counts,one_best)),arb(0)))
        proposal=estimate(upper_arrays(local),K=65536,envelope=UniformInputEnvelope(16,8,4,4),
            occupancies=tuple(range(3,args.max_q+1)),tilt=tilt)
        for q,w in proposal['witnesses'].items():
            m=w['estimated_margin_bits']
            if q not in best or m>best[q]:best[q],choices[q]=m,tilt
        trials.append(proposal)
        result=dict(schema='paired-transvection-screen-1',map_record=record,updates=args.updates,
            sampler='independent u uniform nonzero, v uniform in u-perp including zero; R=I+u*v^T',
            lazy_weight=str(Q(1,2**args.updates)),q1_upper=stride.screen.q1.endpoint(one),
            q1_margin_bits=str(-one.log()/arb(2).log()),q1_fresh_outward=True,precision=192,
            tail_proposal_only=True,whole_code_certificate=False,
            tail_estimated_margin_bits=best,tail_tilt_choices=choices,trials=trials)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        worst=sorted(best,key=best.get)[:5]
        print(f'r={args.updates},tilt={tilt},q1={result["q1_margin_bits"]},worst={[(q,round(best[q],3)) for q in worst]}',flush=True)


if __name__=='__main__':main()
