"""Search fixed quadratic row spaces; this is not setup randomness or a bound.

Use three singleton monomials and six disjoint pairs, covering all fifteen
quadratic monomials on six variables. Edge (0,1) remains a singleton.
Every one of the 512 quadratic combinations is counted exactly by rank.
The best retained map minimizes rank-two then rank-four combinations.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import random

import monomial_maps as maps


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials',type=int,default=20000)
    parser.add_argument('--seed',type=int,default=20261002)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.trials<1 or args.output.exists():
        raise ValueError('positive trial count and fresh output required')
    ranks=[maps.quadratic_rank(mask) for mask in range(1<<15)]
    rng=random.Random(args.seed)
    best=None
    for trial in range(args.trials):
        labels=list(range(1,15)); rng.shuffle(labels)
        basis=[1,1<<labels[0],1<<labels[1]]
        basis += [(1<<labels[i])|(1<<labels[i+1]) for i in range(2,14,2)]
        basis=sorted(basis,key=lambda x:(x.bit_count(),x))
        words=[0]
        for row in basis:
            words += [word^row for word in words]
        counts=Counter(ranks[word] for word in words)
        score=(counts[2],counts[4],tuple(basis))
        if best is None or score<best:
            best=score
            groups=[[list(edge) for i,edge in enumerate(maps.EDGES) if row>>i&1] for row in basis]
            result=dict(proposal_only=True,whole_code_certificate=False,trial=trial,
                trials=args.trials,seed=args.seed,quadratic_basis_masks=basis,
                monomial_groups=groups,quadratic_rank_counts=dict(counts),
                exhaustive_quadratic_census=True,coefficient_xors=0,feedback_xors=6)
            args.output.write_text(json.dumps(result,indent=2)+'\n')
            print(f'trial={trial} rank2={counts[2]} rank4={counts[4]} groups={groups}',flush=True)


if __name__=='__main__':
    main()
