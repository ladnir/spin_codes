"""Bounded fresh screen for alternating uniform-GL and identity updates.

Only M changes: every step emits x+A*a and updates a'=M*a+C*x.
The first physical step in each pair has fresh uniform GL16; the second
has M=I. The identity envelope specializes the retained lazy-refresh
formulas to alpha=1,beta=0. All return rows then use fixed-state fiber
bounds; no fresh-map 1/(2^16-1) return law is applied to an arbitrary state.
The copied updates=0 marker changes no prepared integer map census.
This file never edits or reuses a saved certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import disjoint_pair_maps as maps
import screen_t128 as screen
from flint import arb, arb_mat, ctx
import occupancy_birth_classes as lazy
from packet_outer_geometry_proposal import estimate, upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


def identity_family(physical,z):
    data=dict(physical,updates=0)
    family=lazy.outward_at_z(data,z)
    family=lazy.refine_zero(data,family,z)
    # Capped newborn rows are legitimate here: entering state zero implies
    # a'=Cx for either refresh law. Only that source row is replaced.
    family=screen.q1.kernel_t64.kernel_birth_density.refine_local(data,family,z,Q(1,2))
    n=family[0].nrows(); empty=arb_mat(n,n)
    empty[0,0]=1
    levels=list(map(int,physical['birth_class_levels']))
    empty[1,1]=z**min(levels)
    for i,(level,count) in enumerate(zip(levels,physical['birth_class_counts']),3):
        empty[2,i]=screen.up(screen.aq(Q(int(count),65535))*z**level)
        empty[i,i]=z**level
    family[0]=empty
    return family


def family(data,tilt,identity_first=False):
    physical=screen.q1.kernel_t64.authenticate(data)
    z=(-screen.aq(Q(tilt))).exp()
    fresh=screen.q1.kernel_t64.sparse_kernel.outward_at_z(physical,z)
    fresh=screen.q1.kernel_t64.kernel_birth_density.refine_local(physical,fresh,z,Q(1,2))
    identity=identity_family(physical,z)
    return screen.q1.kernel_t64.convolve(identity,fresh) if identity_first else screen.q1.kernel_t64.convolve(fresh,identity)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--identity-first',action='store_true')
    parser.add_argument('--tilts',nargs='+',required=True)
    parser.add_argument('--max-q',type=int,default=128)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise ValueError('fresh output required')
    ctx.prec=192
    data,record=maps.prepare()
    _,counts=screen.tail.uniform_envelope(16,8,4)
    q1_best=[arb(1)]*65; best={};choices={};trials=[]
    for tilt in args.tilts:
        local=family(data,tilt,args.identity_first)
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=1)
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q1.support_moments(regional[0],regional[1],64)
        q1_best=[min(a,screen.up(factor*b)) for a,b in zip(q1_best,moments)]
        one=screen.up(512*sum((screen.aq(c)*m for c,m in zip(counts,q1_best)),arb(0)))
        proposal=estimate(upper_arrays(local),K=65536,envelope=UniformInputEnvelope(16,8,4,4),
            occupancies=tuple(range(3,args.max_q+1)),tilt=tilt)
        for q,w in proposal['witnesses'].items():
            margin=w['estimated_margin_bits']
            if q not in best or margin>best[q]:best[q],choices[q]=margin,tilt
        trials.append(proposal)
        result=dict(schema='alternating-identity-refresh-screen-1',map_record=record,
            refresh_pattern=['identity','uniform_gl'] if args.identity_first else ['uniform_gl','identity'],
            K=65536,N=131072,threshold=13107,zero_initial_state=True,final_flush=False,
            q1_upper=screen.q1.endpoint(one),q1_margin_bits=str(-one.log()/arb(2).log()),
            q1_fresh_outward=True,precision=192,tail_proposal_only=True,
            whole_code_certificate=False,tail_estimated_margin_bits=best,tail_tilt_choices=choices,trials=trials)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        worst=sorted(best,key=best.get)[:4]
        print(f'tilt={tilt}, q1={result["q1_margin_bits"]}, tailworst={[(q,round(best[q],3)) for q in worst]}',flush=True)


if __name__=='__main__':main()
