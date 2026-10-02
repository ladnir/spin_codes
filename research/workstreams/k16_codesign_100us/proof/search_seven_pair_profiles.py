"""Bounded floating critical-occupancy screen of seven-pair fixed maps.

Every tested candidate gets a fresh full state census and physical kernel.
The critical-q scores are proposals only, not all-occupancy certificates.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path
import random

import monomial_maps as monomial
import seven_pair_maps as fixed
from flint import ctx
from packet_outer_geometry_proposal import estimate,upper_arrays
from rs_uniform_envelope import UniformInputEnvelope

screen=fixed.screen


def prepare(groups):
    rows=[(1<<64)-1]
    rows += [sum(1<<x for x in range(64) if x>>i&1) for i in range(6)]
    for group in groups:
        row=0
        for i,j in group:row ^= sum(1<<x for x in range(64) if (x>>i&1) and (x>>j&1))
        rows.append(row)
    rows=tuple(rows)
    columns=tuple(sum(((row>>x)&1)<<i for i,row in enumerate(rows)) for x in range(64))
    maps=screen.q1.kernel_t64.s16_maps
    if (maps.binary_rank(rows)!=14 or maps.binary_rank(columns)!=14
            or any((a&b).bit_count()&1 for a in rows for b in rows)
            or any(maps.binary_rank(columns[x:x+4])!=4 for x in range(0,64,4))):
        raise ArithmeticError('invalid map')
    images=tuple(maps.images_from_rows(rows));spectrum=Counter(x.bit_count() for x in images)
    if spectrum!={0:1,24:1072,28:3840,32:6558,36:3840,40:1072,64:1}:
        raise ArithmeticError('unexpected spectrum')
    physical=screen.q1.kernel_t64.kernel_maps.prepare_maps(images,columns,bits=14,
        distribution='uniform_gl',birth_density='capped')
    return screen.q1.kernel_t64.wrap(physical),dict(monomial_groups=groups,
        expansion_rows_hex=list(map(hex,rows)),feedback_columns=columns,
        spectrum={str(w):n for w,n in sorted(spectrum.items())},
        map_sha256=physical['map_sha256'],full_state_census=True)


def candidates(count,seed):
    yield fixed.GROUPS
    ranks=[monomial.quadratic_rank(mask) for mask in range(1<<15)]
    compatible={i:tuple(j for j in range(15) if not set(monomial.EDGES[i])&set(monomial.EDGES[j])) for i in range(15)}
    rng=random.Random(seed);seen={fixed.GROUPS};accepted=1
    while accepted<count:
        omitted=rng.randrange(1,15)
        remaining=list(range(15));remaining.remove(omitted);basis=[]
        while remaining:
            i=remaining.pop(rng.randrange(len(remaining)))
            choices=[j for j in remaining if j in compatible[i]]
            if not choices:break
            j=rng.choice(choices);remaining.remove(j);basis.append((1<<i)|(1<<j))
        if remaining or len(basis)!=7:continue
        basis.sort();words=[0]
        for row in basis:words += [word^row for word in words]
        if any(ranks[word]==2 for word in words):continue
        groups=tuple(tuple(edge for i,edge in enumerate(monomial.EDGES) if row>>i&1) for row in basis)
        if groups in seen:continue
        seen.add(groups);accepted+=1;yield groups


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates',type=int,default=16)
    parser.add_argument('--seed',type=int,default=20261004)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.candidates<1 or args.output.exists():raise ValueError('positive count and fresh output required')
    ctx.prec=160;trials=[]
    for index,groups in enumerate(candidates(args.candidates,args.seed)):
        data,record=prepare(groups);best={};choices={}
        for tilt in ('.13','.1325','.135','.14','.1425','.145','.1475','.15','.155','.16','.1625','.165'):
            local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
            result=estimate(upper_arrays(local),K=65536,envelope=UniformInputEnvelope(16,8,4,4),
                occupancies=(36,38,40,42,44),tilt=tilt)
            for q,w in result['witnesses'].items():
                margin=w['estimated_margin_bits']
                if q not in best or margin>best[q]:best[q],choices[q]=margin,tilt
        score=min(best.values())
        trials.append(dict(index=index,map_record=record,critical_margin_bits=best,
            critical_tilt_choices=choices,score=score))
        ordered=sorted(trials,key=lambda r:r['score'],reverse=True)
        args.output.write_text(json.dumps(dict(schema='seven-pair-critical-profile-search-1',
            proposal_only=True,whole_code_certificate=False,seed=args.seed,
            requested_candidates=args.candidates,trials=trials,best=ordered[0]),indent=2)+'\n')
        print(f'candidate={index}, score={score:.6f}, margins={best}, best={ordered[0]["index"]}',flush=True)


if __name__=='__main__':main()
