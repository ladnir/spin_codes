"""Fresh partial q1 and floating tail screen for the fixed seven-pair map."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import seven_pair_maps as maps
from flint import arb,ctx
from packet_outer_geometry_proposal import estimate,upper_arrays
from rs_uniform_envelope import UniformInputEnvelope

screen=maps.screen


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',required=True)
    parser.add_argument('--max-q',type=int,default=128)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('fresh output required')
    ctx.prec=192;data,record=maps.prepare()
    _,counts=screen.tail.uniform_envelope(16,8,4)
    one_best=[arb(1)]*65;best={};choices={};trials=[]
    for tilt in args.tilts:
        local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=1)
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q1.support_moments(regional[0],regional[1],64)
        one_best=[min(a,screen.up(factor*b)) for a,b in zip(one_best,moments)]
        one=screen.up(512*sum((screen.aq(c)*m for c,m in zip(counts,one_best)),arb(0)))
        proposal=estimate(upper_arrays(local),K=65536,envelope=UniformInputEnvelope(16,8,4,4),
            occupancies=tuple(range(3,args.max_q+1)),tilt=tilt)
        for q,w in proposal['witnesses'].items():
            margin=w['estimated_margin_bits']
            if q not in best or margin>best[q]:best[q],choices[q]=margin,tilt
        trials.append(proposal)
        result=dict(schema='seven-pair-map-screen-1',map_record=record,
            K=65536,N=131072,threshold=13107,zero_initial_state=True,final_flush=False,
            q1_upper=screen.q1.endpoint(one),q1_margin_bits=str(-one.log()/arb(2).log()),
            q1_fresh_outward=True,precision=192,tail_proposal_only=True,
            whole_code_certificate=False,tail_estimated_margin_bits=best,tail_tilt_choices=choices,trials=trials)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        worst=sorted(best,key=best.get)[:4]
        print(f'tilt={tilt}, q1={result["q1_margin_bits"]}, worst={[(q,round(best[q],3)) for q in worst]}',flush=True)


if __name__=='__main__':main()
