"""Bounded proposal search for s14: seven affine and seven paired rows.

The quadratic rows partition fourteen of the fifteen quadratic monomials.
Each pair uses vertex-disjoint edges, so every individual row has rank four.
The omitted edge is never (0,1), preserving rank-four packet restrictions.
This exact rank census selects a candidate; it is not a distance certificate.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import random

import monomial_maps as maps


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials',type=int,default=30000)
    parser.add_argument('--seed',type=int,default=20261003)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.trials<1 or args.output.exists():raise ValueError('positive trials and fresh output required')
    ranks=[maps.quadratic_rank(mask) for mask in range(1<<15)]
    compatible={i:tuple(j for j in range(15) if not set(maps.EDGES[i])&set(maps.EDGES[j])) for i in range(15)}
    rng=random.Random(args.seed);best=None;accepted=0;attempts=0
    while accepted<args.trials:
        attempts+=1
        omitted=rng.randrange(1,15)
        remaining=list(range(15));remaining.remove(omitted);basis=[]
        while remaining:
            i=remaining.pop(rng.randrange(len(remaining)))
            choices=[j for j in remaining if j in compatible[i]]
            if not choices:break
            j=rng.choice(choices);remaining.remove(j);basis.append((1<<i)|(1<<j))
        if remaining or len(basis)!=7:continue
        accepted+=1;basis.sort();words=[0]
        for row in basis:words += [word^row for word in words]
        counts=Counter(ranks[word] for word in words)
        score=(counts[2],counts[4],tuple(basis))
        if best is None or score<best:
            best=score
            groups=[[list(edge) for i,edge in enumerate(maps.EDGES) if row>>i&1] for row in basis]
            result=dict(schema='s14-seven-pair-proposal-1',proposal_only=True,
                whole_code_certificate=False,accepted_trial=accepted,attempts=attempts,
                trials=args.trials,seed=args.seed,quadratic_basis_masks=basis,
                omitted_edge=list(maps.EDGES[omitted]),monomial_groups=groups,
                quadratic_rank_counts=dict(counts),exhaustive_quadratic_census=True,
                coefficient_xors=0,feedback_xors=7)
            args.output.write_text(json.dumps(result,indent=2)+'\n')
            print(f'trial={accepted}, rank2={counts[2]}, rank4={counts[4]}, groups={groups}',flush=True)


if __name__=='__main__':main()
