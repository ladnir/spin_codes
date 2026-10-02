"""Floating proposals for one identity refresh in each fixed period.

Chronological occupancy convolutions retain every physical step and its
input distribution. No state reset occurs at period or region boundaries.
Local matrices come from valid fresh-GL and identity envelopes, but their
floating compositions are proposals and are not outward bounds.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path

import numpy as np
from flint import ctx
import screen_refresh_stride as stride
from packet_outer_geometry_proposal import estimate,upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


def convolve(left,right):
    W,V=len(left)-1,len(right)-1
    result=np.zeros((W+V+1,left.shape[1],left.shape[2]))
    for j in range(W+V+1):
        denominator=comb(W+V,j)
        for a in range(max(0,j-V),min(W,j)+1):
            weight=comb(W,a)*comb(V,j-a)/denominator
            result[j]+=weight*(left[a]@right[j-a])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--period',type=int,choices=(4,8,16),required=True)
    parser.add_argument('--tilts',nargs='+',required=True)
    parser.add_argument('--max-q',type=int,default=128)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('fresh output required')
    ctx.prec=192
    wrapped,record=stride.maps.prepare()
    physical=stride.screen.q1.kernel_t64.authenticate(wrapped)
    best={};choices={};trials=[]
    for tilt in args.tilts:
        z=(-stride.screen.aq(Q(tilt))).exp()
        fresh=stride.screen.q1.kernel_t64.sparse_kernel.outward_at_z(physical,z)
        fresh=stride.screen.q1.kernel_t64.kernel_birth_density.refine_local(physical,fresh,z,Q(1,2))
        fresh=upper_arrays(fresh)
        identity=upper_arrays(stride.identity_family(physical,z))
        local=fresh
        with np.errstate(over='raise',invalid='raise',under='raise'):
            for step in range(1,args.period-1):local=convolve(local,fresh)
            local=convolve(local,identity)
        proposal=estimate(local,K=65536,envelope=UniformInputEnvelope(16,8,4,4),
            occupancies=tuple(range(3,args.max_q+1)),tilt=tilt,windows=16*args.period)
        for q,w in proposal['witnesses'].items():
            m=w['estimated_margin_bits']
            if q not in best or m>best[q]:best[q],choices[q]=m,tilt
        trials.append(proposal)
        args.output.write_text(json.dumps(dict(proposal_only=True,whole_code_certificate=False,
            period=args.period,refresh_pattern=['uniform_gl']*(args.period-1)+['identity'],
            map_record=record,tail_estimated_margin_bits=best,tail_tilt_choices=choices,trials=trials),indent=2)+'\n')
        worst=sorted(best,key=best.get)[:5]
        print(f'period={args.period},tilt={tilt},worst={[(q,round(best[q],3)) for q in worst]}',flush=True)


if __name__=='__main__':main()
