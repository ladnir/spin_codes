"""Adaptive complete support covers for the state-memory envelope.

Boxes are accepted only against a volume-proportional failure budget.
Splitting repeated intervals retains exact label multiplicities. Numerical
witness selection is separate from optional outward replay.
"""
import argparse
from collections import Counter
from math import comb,log,prod
import numpy as np
from scipy.special import logsumexp
from flint import arb,ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from group_rank_one_verify import up
from occupancy_model import local_data
from occupancy_memory import prepare,build,TERMINAL
from occupancy_memory_verify import TILTS,DENOMINATOR,replay
from occupancy_screen import optimize,matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass


def volume(intervals,multiplicity):
    return multiplicity*prod(hi-lo+1 for lo,hi in intervals)


def split(intervals,multiplicity):
    chosen=max(intervals,key=lambda p:(p[1]-p[0],-p[0]))
    lo,hi=chosen
    assert lo<hi
    middle=(lo+hi)//2
    repeats=intervals.count(chosen)
    other=[p for p in intervals if p!=chosen]
    children=[]
    for k in range(repeats+1):
        child=tuple(sorted(other+[(lo,middle)]*k+[(middle+1,hi)]*(repeats-k)))
        children.append((child,multiplicity*comb(repeats,k)))
    assert sum(volume(*c) for c in children)==volume(intervals,multiplicity)
    return children


def geometry_test():
    from itertools import product,permutations
    for q in range(1,5):
        pending=[(((0,4),)*q,1)]
        covered=Counter()
        while pending:
            box,mult=pending.pop()
            if any(hi-lo>1 for lo,hi in box):
                pending.extend(split(box,mult))
                continue
            # At this partition each box's permutations are disjoint and
            # the stored multiplicity equals the number of labelings.
            labelings=set(permutations(box))
            assert len(labelings)==mult
            for labeled in labelings:
                covered.update(product(*(range(lo,hi+1) for lo,hi in labeled)))
        assert len(covered)==5**q and set(covered.values())=={1}
    print('Adaptive box geometry and exact label multiplicities checked',flush=True)


def witness_selector(operators,groups,balanced=False):
    cache={}
    diagonal_cache={}
    def diagonal(tilt,u):
        key=tilt,u
        if key not in diagonal_cache:
            _,ps=optimize(operators[tilt][1],(u,)*groups,float(tilt),TERMINAL)
            diagonal_cache[key]=ps[0]
        return diagonal_cache[key]
    def probability(tilt,interval,policy):
        key=tilt,interval,policy
        if key not in cache:
            lo,hi=interval
            p=diagonal(tilt,(lo+hi)/2)
            if policy and lo<hi:
                # Equal endpoint binomial masses, clamped between the two
                # diagonal point witnesses. This only proposes a witness;
                # validity never depends on its optimality.
                z=(log(comb(256,lo))-log(comb(256,hi)))/(hi-lo)
                equal=1/(1+np.exp(-z))
                low,high=sorted((diagonal(tilt,lo),diagonal(tilt,hi)))
                p=max(low,min(high,equal))
            cache[key]=DENOMINATOR if lo==256 else max(1,min(DENOMINATOR-1,round(p*DENOMINATOR)))
        return cache[key]
    def witness(box):
        best=None
        for tilt,(_,region) in operators.items():
            for policy in ((False,True) if balanced else (False,)):
                nums=[probability(tilt,p,policy) for p in box]
                ps=[n/DENOMINATOR for n in nums]
                value=log_power_moment(matrix_for_probabilities(region,ps),256,TERMINAL)+float(tilt)*209715
                value-=sum(min(log_binomial_mass(256,u,p) for u in (lo,hi)) for p,(lo,hi) in zip(ps,box))
                if best is None or value<best[0]:
                    best=value,tilt,nums
        return best
    return witness


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=5,choices=range(3,33))
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--target-bits',type=int,default=55)
    parser.add_argument('--max-nodes',type=int,default=50000)
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--tilts',nargs='+',default=list(TILTS)+['.0032','.004','.005','.0064','.008','.01'])
    args=parser.parse_args()
    assert args.precision>=128
    geometry_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    counts=[sum(row[u] for row in caps) for u in range(257)]
    prepared=prepare(local_data(min(args.groups,4)))
    ctx.prec=args.precision
    operators={}
    for tilt in args.tilts:
        operators[tilt]=build(prepared,tilt,args.groups)
        print('Built adaptive operators',tilt,flush=True)
    witness=witness_selector(operators,args.groups)
    total_volume=219**args.groups
    pending=[(((38,256),)*args.groups,1)]
    accepted=0
    accepted_volume=0
    terms=[]
    total=arb(0)
    visited=0
    unresolved=[]
    location_count=comb(2048,args.groups)
    while pending and visited<args.max_nodes:
        box,mult=pending.pop()
        visited+=1
        size=volume(box,mult)
        count=location_count*mult*prod(counts[hi] for lo,hi in box)
        value,tilt,ps=witness(box)
        score=min(0.,value)+log(count)
        budget=-args.target_bits*log(2)+log(size)-log(total_volume)
        if score<budget-1e-6:
            if not args.screen_only:
                rectangles=[(lo,hi,counts[hi]) for lo,hi in box]
                term=up(replay(operators[tilt][0],tilt,ps,rectangles)*count)
                allowed=arb(2)**(-args.target_bits)*size/total_volume
                assert term<allowed,('outward acceptance failed',box)
                total=up(total+term)
            accepted+=1
            accepted_volume+=size
            terms.append(score)
        elif all(lo==hi for lo,hi in box):
            unresolved.append((score/log(2),box,mult))
        else:
            pending.extend(split(box,mult))
        if visited%1000==0:
            print('Adaptive nodes',visited,'accepted',accepted,'pending',len(pending),'unresolved',len(unresolved),flush=True)
    missing=sum(volume(box,mult) for box,mult in pending)+sum(volume(box,mult) for _,box,mult in unresolved)
    assert accepted_volume+missing==total_volume
    print('Adaptive coverage:',visited,'nodes;',accepted,'accepted leaves;',len(pending),'pending;',len(unresolved),'unresolved',flush=True)
    print('BINARY64 accepted union log2',logsumexp(terms)/log(2) if terms else '-inf',flush=True)
    if unresolved:
        print('Largest unresolved points:',*sorted(unresolved,reverse=True)[:5],sep='\n',flush=True)
    if missing:
        print('INCOMPLETE: missing support volume',missing,'of',total_volume,flush=True)
    elif not args.screen_only:
        assert 0<total<arb(2)**-args.target_bits
        print('VERIFIED occupancy',args.groups,'precision',args.precision,'upper',total,
              'margin',-total.log()/arb(2).log(),flush=True)
    print('Full-code certificate: NO. Other occupancies remain.')


if __name__=='__main__':
    main()
