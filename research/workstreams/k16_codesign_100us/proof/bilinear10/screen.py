"""Fresh partial q1 and critical-tail screen; never a whole-code certificate."""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
from time import monotonic

import bilinear_maps as maps
from flint import arb,ctx
from packet_outer_geometry_proposal import estimate,upper_arrays
from rs_uniform_envelope import UniformInputEnvelope

screen=maps.screen


def run(tilts,tail_q,outward_q,output,precision):
    output=Path(output).resolve()
    if output.exists():raise FileExistsError('fresh output required')
    if precision<192 or not tilts or any(Q(t)<=0 for t in tilts):raise ValueError('positive tilts and precision>=192 required')
    tail_q=tuple(sorted(set(tail_q)));outward_q=tuple(sorted(set(outward_q)))
    if any(not 3<=q<=512 for q in tail_q+outward_q):raise ValueError('tail occupancy must be3..512')
    ctx.prec=precision;start=monotonic()
    data,map_record=maps.prepare();pins=maps.source_pins()
    beta,counts=screen.tail.uniform_envelope(16,8,4)
    one=[arb(1)]*65;best={};choices={};exact={};exact_choices={};trials=[]
    result=dict(schema='bilinear10-partial-screen-1',K=65536,N=131072,threshold=13107,
        outer='four GF16 RS[16,8] rows with independent transitive16-bit symbol maps',
        groups=512,group_dimension=128,regions=64,packet_bits=4,
        physical_steps_per_region=32,macro_steps_per_region=16,
        zero_initial_state=True,final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        routing='independent uniform within-group packet and regional slot permutations',
        map_record=map_record,source_sha256=pins,precision=precision,
        tilts=list(tilts),q1_fresh_outward=True,tail_proposal_only=True,
        whole_code_certificate=False,
        scope='Fresh partial bounds for the explicit ideal independent GL10 ensemble. No full certificate, deterministic-seed guarantee, or timing claim.')
    print('Fresh algebra:',json.dumps(map_record['expansion_spectrum']),flush=True)
    for tilt in tilts:
        local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,
                                     maximum_groups=max((1,)+outward_q))
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q1.support_moments(regional[0],regional[1],64)
        one=[min(a,screen.up(factor*b)) for a,b in zip(one,moments)]
        q1=screen.up(512*sum((screen.aq(c)*m for c,m in zip(counts,one)),arb(0)))
        if tail_q:
            proposal=estimate(upper_arrays(local),K=65536,envelope=UniformInputEnvelope(16,8,4,4),
                              occupancies=tail_q,tilt=tilt)
            trials.append(proposal)
            for q,w in proposal['witnesses'].items():
                margin=w['estimated_margin_bits']
                if q not in best or margin>best[q]:best[q],choices[q]=margin,tilt
        for q in outward_q:
            matrix=screen.tail.regional_uniform(regional,q)**64
            moment=sum((matrix[0,j] for j in range(matrix.ncols())),arb(0))
            value=screen.up(comb(512,q)*screen.aq(beta)**q*factor*moment)
            if q not in exact or value<exact[q]:exact[q],exact_choices[str(q)]=value,tilt
        if pins!=maps.source_pins() or ctx.prec!=precision:raise ArithmeticError('source or precision changed')
        result.update(q1_upper=screen.q1.endpoint(q1),q1_margin_bits=str(-q1.log()/arb(2).log()),
            q1_support_probability_uppers=[screen.q1.endpoint(x) for x in one],
            tail_estimated_margin_bits=best,tail_tilt_choices=choices,tail_trials=trials,
            outward_occupancy_uppers={str(q):screen.q1.endpoint(x) for q,x in exact.items()},
            outward_occupancy_margin_bits={str(q):str(-x.log()/arb(2).log()) for q,x in exact.items()},
            outward_choices=exact_choices,elapsed_seconds=monotonic()-start)
        output.write_text(json.dumps(result,indent=2)+'\n')
        worst=sorted(best,key=best.get)[:3]
        print('tilt',tilt,'q1',result['q1_margin_bits'],'tailworst',[(q,round(best[q],5)) for q in worst],flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.00064','.00128','.00256','.00512','.01024','.0256','.0512','.1024','.2048','.4096','.8192'])
    parser.add_argument('--tail-q',type=int,nargs='*',default=[3,8,16,24,32,48,64,80,96,112,128])
    parser.add_argument('--outward-q',type=int,nargs='*',default=[])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.tilts,args.tail_q,args.outward_q,args.output,args.precision)
