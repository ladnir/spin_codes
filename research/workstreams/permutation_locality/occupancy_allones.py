"""All-one column counts coupled to the nine-coordinate inner envelope.

Selected-point binary64 diagnostic only, not a full-code certificate.
"""
import argparse
from functools import lru_cache
from itertools import product
from math import comb,log,prod
import numpy as np
from flint import ctx

from bch_joint_support import authenticated_caps,support_caps,rank_total
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_count_gap import oa_refinement
from occupancy_weighted import log_weighted_cap
from occupancy_memory import prepare,epoch_operators,rounded,TERMINAL
from occupancy_model import local_data,placement
from occupancy_screen import optimize
from joint_support import span,gf2_rank


def full_column_caps(spectrum,caps,dimensions,u,total_cap):
    """Upper counts by number j of all-one columns, union support exactly u.

    For j>0, even combinations span a codimension-one subspace H.
    Its support is v=u-j. Each extension coset gives affine spanning
    four-tuples counted by 2^(h-1) prod_{i<h-1}(8-2^i).
    """
    n=len(spectrum)-1
    result=[total_cap]+[0]*u
    cdf=[sum(spectrum[1:w+1]) for w in range(n+1)]
    for j in range(1,u+1):
        v=u-j
        # Rank one: the sole tuple with an all-one column is (w,w,w,w).
        value=spectrum[u] if v==0 else 0
        for h in range(2,5):
            if dimensions[v]<h-1:
                continue
            spaces=caps[h-2][v]//rank_total(h-1,4,h-1)
            extensions=min((comb(n-v,j)*(1<<dimensions[v]))//(1<<(h-1)),cdf[(2*u-v)//2])
            tuples=(1<<(h-1))*prod(8-(1<<i) for i in range(h-1))
            value+=min(caps[h-1][u],tuples*spaces*extensions)
        result[j]=min(total_cap,value)
    return result


@lru_cache(maxsize=4)
def prefix_flag_caps(spectrum,caps,dimensions):
    """Counts by all-one count J and union support <=u, using H's CDF once."""
    n=len(spectrum)-1
    cdf=[sum(spectrum[1:w+1]) for w in range(n+1)]
    result=[[0]*(n+1) for _ in range(n+1)]
    for j in range(1,n+1):
        for u in range(j,n+1):
            result[u][j]=spectrum[j]  # rank-one tuple (w,w,w,w)
        for h in range(2,5):
            denominator=rank_total(h-1,4,h-1)
            spaces=[row//denominator for row in caps[h-2]]
            increments=[spaces[0]]+[spaces[v]-spaces[v-1] for v in range(1,n+1)]
            tuples=(1<<(h-1))*prod(8-(1<<i) for i in range(h-1))
            extensions=[min((comb(n-v,j)*(1<<dimensions[v]))//(1<<(h-1)),cdf[(v+2*j)//2])
                        if dimensions[v]>=h-1 else 0 for v in range(n-j+1)]
            for u in range(j,n+1):
                maximum=0;bound=0
                for v in range(u-j,-1,-1):
                    maximum=max(maximum,extensions[v])
                    bound+=increments[v]*maximum
                result[u][j]+=min(caps[h-1][u],tuples*bound)
    return tuple(tuple(row) for row in result)


def weighted_cdf(spectrum,caps,dimensions,penalty,prefix_flags=False,extra_prefix=None):
    """Integer upper CDF for counts weighted by penalty^-J."""
    from fractions import Fraction
    from occupancy_weighted import allocate
    rho=Fraction(penalty)
    assert 0<rho<=1
    totals=[sum(row[u] for row in caps) for u in range(len(spectrum))]
    if rho==1:
        return totals
    prefix=(prefix_flag_caps(tuple(spectrum),tuple(tuple(row) for row in caps),tuple(dimensions))
            if prefix_flags else None)
    shells=[0]*len(spectrum)
    result=[0]
    for u in range(1,len(spectrum)):
        current=full_column_caps(spectrum,caps,dimensions,u,totals[u])
        for j,value in enumerate(current):
            shells[j]+=value
        bounds=([min(value,prefix[u][j]) if j else value for j,value in enumerate(shells)]
                if prefix is not None else shells)
        if extra_prefix is not None:
            bounds=[min(value,extra_prefix[u][j]) if j else value for j,value in enumerate(bounds)]
        upper=sum((c*rho**(-j) for j,c in allocate(bounds,totals[u],rho)),Fraction(0))
        result.append(-(-upper.numerator//upper.denominator))
    assert all(a<=b for a,b in zip(result,result[1:]))
    return result


def self_test():
    checks=0
    for rows in ([1,2,4],[0b1111000,0b1100110,0b1010101],[0x97,0x4b,0x2d,0x1e]):
        words=span(rows)
        n=max(words).bit_length()
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        ranks=[[0]*(n+1) for _ in range(4)]
        actual=[[0]*(u+1) for u in range(n+1)]
        for tup in product(words,repeat=4):
            u=(tup[0]|tup[1]|tup[2]|tup[3]).bit_count()
            j=(tup[0]&tup[1]&tup[2]&tup[3]).bit_count()
            actual[u][j]+=1
            h=gf2_rank(tup)
            if h:
                ranks[h-1][u]+=1
        caps=[]
        for row in ranks:
            caps.append([sum(row[:u+1]) for u in range(n+1)])
        dimensions=[0]*(n+1)
        for mask in range(1<<n):
            size=sum(w&~mask==0 for w in words)
            dimensions[mask.bit_count()]=max(dimensions[mask.bit_count()],size.bit_length()-1)
        for u in range(1,n+1):
            upper=full_column_caps(spectrum,caps,dimensions,u,sum(actual[u]))
            assert all(a<=b for a,b in zip(actual[u],upper))
            checks+=len(upper)
        from fractions import Fraction
        for rho in (Fraction(1,2),Fraction(3,4),Fraction(1)):
            upper=weighted_cdf(spectrum,caps,dimensions,rho)
            improved=weighted_cdf(spectrum,caps,dimensions,rho,prefix_flags=True)
            cumulative=Fraction(0)
            for u in range(1,n+1):
                cumulative+=sum(c*rho**(-j) for j,c in enumerate(actual[u]))
                assert cumulative<=upper[u]
                assert cumulative<=improved[u]<=upper[u]
    print('All-one-column flag counts:',checks,'exact small-code shell inequalities',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=16,choices=range(4,33))
    parser.add_argument('--fresh',action='store_true')
    parser.add_argument('--pair-memory',action='store_true')
    parser.add_argument('--window-average',action='store_true')
    parser.add_argument('--multi-average',action='store_true')
    parser.add_argument('--tilts',nargs='+',default=['.005','.0064','.008','.01','.0128'])
    parser.add_argument('--penalties',nargs='+',default=['.5','.75','.9','1','1.1'])
    parser.add_argument('--supports',nargs='+',type=int,default=[76,80,84,96,128,176,216])
    args=parser.parse_args()
    assert not(args.fresh and args.pair_memory)
    assert not(args.window_average and args.pair_memory)
    assert not(args.multi_average and args.pair_memory)
    self_test()
    spectrum=authenticated_caps()
    dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    oa,_=oa_refinement(caps)
    totals={u:sum(row[u] for row in oa) for u in args.supports}
    shells={u:full_column_caps(spectrum,caps,dimensions,u,totals[u]) for u in args.supports}
    counts={(u,p):log_weighted_cap(shells[u],totals[u],float(p)) for u in args.supports for p in args.penalties}
    prepared=prepare(local_data(4))
    if args.fresh or args.multi_average:
        from occupancy_fresh_moment import fresh_census,refine
        census=fresh_census(prepared)
    if args.multi_average:
        import occupancy_multi_average as multi_average
        multi_average.self_test()
    terminal=TERMINAL
    if args.pair_memory:
        import occupancy_pair_memory as pair_memory
        census=pair_memory.census(prepared)
        terminal=pair_memory.PAIR_TERMINAL
    if args.window_average:
        import occupancy_window_average as window_average
        windows=window_average.prepare_inputs()
    ctx.prec=192
    if args.pair_memory:
        pair_memory.self_test(prepared,census)
    best={u:(float('inf'),None,None) for u in args.supports}
    for tilt in args.tilts:
        if args.window_average:
            averages=window_average.averages(windows,tilt)
            window_average.self_test(windows,averages,tilt)
        for penalty in args.penalties:
            if args.pair_memory:
                ops=pair_memory.operators(prepared,census,tilt,args.groups,full_penalty=penalty)
            elif args.fresh:
                ops=refine(prepared,census,tilt,args.groups,full_penalty=penalty)
            else:
                ops=epoch_operators(prepared,tilt,maximum_groups=args.groups,pair_conditioned=True,
                                    full_penalty=penalty if penalty!='1' else 1)
            if args.window_average:
                ops=window_average.refine(prepared,averages,tilt,args.groups,penalty,ops)
            if args.multi_average:
                ops=multi_average.refine(ops,census,tilt,penalty)
            regions=placement(ops,rounding=rounded)
            size=len(terminal)
            arrays=[np.array([[float(t[i,j]) for j in range(size)] for i in range(size)]) for t in regions]
            for u in args.supports:
                value,_=optimize(arrays,(u,)*args.groups,float(tilt),terminal)
                score=(value+args.groups*counts[u,penalty]+log(comb(2048,args.groups)))/log(2)
                if score<best[u][0]:
                    best[u]=float(score),tilt,penalty
            print('BINARY64 all-one-column screen',tilt,penalty,'best',[(u,round(best[u][0],3)) for u in args.supports],flush=True)
    for u,row in best.items():
        print('RESULT u, log2 contribution, tilt, penalty',u,*row,flush=True)
    print('Selected supports only; OA premise required; no full-code certificate.')


if __name__=='__main__':
    main()
