"""Fresh full K16 replay after deleting paired state row10; GL15 updates.

The outer/routing geometry is the completed paired-map construction.
This candidate removes x0*x2+x1*x4 from BOTH expansion and transpose
feedback, and samples a new independent uniform GL15 map every t64 step.
No prior numerical endpoint is reused. Source modules of the s16 certificate
remain unchanged. The rational tilt recipe is reused only as a search hint.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
from time import monotonic

from flint import arb,ctx
import screen_paired_restriction as restriction
import reproduce_disjoint_pairs as shared

screen=shared.screen


def run(output):
    output=Path(output).resolve()
    if output.exists():raise FileExistsError('fresh output required')
    ctx.prec=256;start=monotonic()
    data,map_record=restriction.prepare((10,))
    beta,counts=screen.tail.uniform_envelope(16,8,4)
    pins=shared.sources();values={};choices={}
    record=dict(schema='k16-paired-s15-drop10-complete-replay-1',K=65536,N=131072,
        threshold=13107,distance='1/10',target_margin_bits=40,target_minimum_distance=13108,
        groups=512,group_dimension=128,regions=64,packet_bits=4,
        physical_t=64,state_bits=15,physical_steps=2048,physical_steps_per_region=32,
        zero_initial_state=True,final_flush=False,state_continuity='retained_across_every_step_and_region',
        outer='four parallel GF16 RS[16,8] rows with independent transitive GL16 symbol maps',
        inner_randomness='independent uniform GL15 matrix for each physical t64 step',
        routing_randomness='independent uniform per-group packet and regional slot permutations',
        map_record=map_record,source_sha256=pins,precision=256,beta=str(beta),
        q1_tilts=shared.Q1_TILTS,q2_tilts=shared.Q2_TILTS,tail_recipes=shared.TAIL_RECIPES,
        fresh_replay=True,all_occupancies_covered=False,whole_code_certificate=False,
        scope='First-moment setup-failure bound for the explicit ideal GL15 ensemble; no GF(2^16) state-family transfer or deterministic-seed guarantee.')

    def checkpoint():
        if pins!=shared.sources() or ctx.prec!=256:
            raise ArithmeticError('source identity or precision changed')
        if any(not x.is_finite() or not x>0 for x in values.values()):
            raise ArithmeticError('positive finite endpoints required')
        record.update(occupancy_uppers={str(q):screen.q1.endpoint(x) for q,x in values.items()},
            occupancy_margin_bits={str(q):shared.margin(x) for q,x in values.items()},
            tail_choices=choices,elapsed_seconds=monotonic()-start)
        output.write_text(json.dumps(record,indent=2)+'\n')

    one=[arb(1)]*65
    for tilt in shared.Q1_TILTS:
        local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=1)
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q1.support_moments(regional[0],regional[1],64)
        one=[min(a,screen.up(factor*b)) for a,b in zip(one,moments)]
    values[1]=screen.up(512*sum((screen.aq(c)*m for c,m in zip(counts,one)),arb(0)))
    record['q1_support_probability_uppers']=[screen.q1.endpoint(x) for x in one]
    checkpoint();print(f'q1={shared.margin(values[1])}',flush=True)
    two=[[arb(1)]*(v+1) for v in range(65)]
    for tilt in shared.Q2_TILTS:
        local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=2)
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q2.pair_support_moments(regional,64)
        for v,column in enumerate(moments):
            for u,m in enumerate(column):two[v][u]=min(two[v][u],screen.up(factor*m))
    values[2]=screen.q2.fold_shell_pairs(counts,two,512)
    checkpoint();print(f'q2={shared.margin(values[2])}',flush=True)
    for first,last,tilts in shared.TAIL_RECIPES:
        for tilt in tilts:
            local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
            regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=last)
            factor=(screen.aq(Q(tilt))*13107).exp()
            for q in range(first,last+1):
                matrix=screen.tail.regional_uniform(regional,q)**64
                moment=sum((matrix[0,j] for j in range(matrix.ncols())),arb(0))
                value=screen.up(comb(512,q)*screen.aq(beta)**q*factor*moment)
                if q not in values or value<values[q]:values[q],choices[str(q)]=value,tilt
            checkpoint()
            worst=min(range(first,last+1),key=lambda q:float(-values[q].log()/arb(2).log()))
            print(f'q={first}..{last}, tilt={tilt}, worst q{worst}={shared.margin(values[worst])}',flush=True)
    if set(values)!=set(range(1,513)):raise ArithmeticError('incomplete coverage')
    total=shared.exact_sum.dyadic_sum(screen.q1.endpoint(v) for v in values.values())
    endpoint=shared.exact_sum.rounded_endpoint(total,256)
    passed=shared.exact_sum.dyadic_less(shared.exact_sum.dyadic(endpoint),(1,-40))
    record.update(all_occupancies_covered=True,occupancy_covered=[1,512],target_met=passed,
        whole_code_certificate=passed,union_exact_mantissa_hex=hex(total[0]),union_exact_exponent=total[1],
        union_upper=endpoint,margin_bits=shared.margin(arb(endpoint[0])*arb(2)**endpoint[1]),
        arithmetic='Exact dyadic endpoint sum, one upward rounding, exact integer comparison to target.')
    checkpoint();print(f'Complete union={record["margin_bits"]}, target_met={passed}',flush=True)
    if not passed:raise ArithmeticError('fresh union did not meet target')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
