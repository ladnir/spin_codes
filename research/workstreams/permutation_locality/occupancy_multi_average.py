"""Without-replacement multi-window moment bound, preserving the nine-state cone.

Only mass/refresh entries and an explicit uniform-class zero upper change.
Density entries are not changed by this independence argument.
"""
from fractions import Fraction
from itertools import combinations_with_replacement,permutations,product
from math import comb,prod
from flint import arb

from group_rank_one_verify import up
from occupancy_memory import Z,F,M,C,U


def multiplicities(j):
    for n1 in range(j+1):
        for n2 in range(j-n1+1):
            for n3 in range(j-n1-n2+1):
                yield (n1,n2,n3,j-n1-n2-n3)


def self_test():
    checks=0
    width=4
    for z in (Fraction(7,8),Fraction(15,16)):
        phi={(a,b):sum((Fraction(comb(b,k)*comb(width-b,a-k),comb(width,a))*z**(-2*k)
                           for k in range(max(0,a+b-width),min(a,b)+1)),Fraction(0))
             for a in range(1,5) for b in range(5)}
        assert all(phi[a,b]<=phi[a,b+1] for a in range(1,5) for b in range(4))
        for n in (2,3,4):
            for population in combinations_with_replacement(range(5),n):
                v=sum(population)
                mean={a:sum(phi[a,b] for b in population)/n for a in range(1,5)}
                for j in range(1,n+1):
                    choices=list(permutations(range(n),j))
                    # Sorted weights suffice: assignment is uniform over
                    # labeled distinct windows, so the expectation is symmetric.
                    for weights in combinations_with_replacement(range(1,5),j):
                        actual=sum((prod(phi[a,population[k]] for a,k in zip(weights,choice))
                                    for choice in choices),Fraction(0))/len(choices)
                        independent=prod(mean[a] for a in weights)
                        chord=prod(1+Fraction(v,n*width)*(z**(-2*a)-1) for a in weights)
                        assert actual<=independent<=chord
                        checks+=1
    print('Exact distinct-window product and chord checks:',checks,flush=True)


def refine(base,fresh_census,tilt,full_penalty=1,input_penalty=1,odd_penalty=1):
    spectrum,distributions,_=fresh_census
    levels=sorted(spectrum)
    m=(1<<19)-1
    rho=arb(full_penalty)
    factors={v:[(arb(128-v)*(-arb(tilt)*a).exp()+arb(v)*(arb(tilt)*a).exp())/128
                for a in range(1,5)] for v in levels}
    initial={v:(-arb(tilt)*v).exp() for v in levels}
    affected=[(F,M),(M,M),(M,Z)]
    affected += [(i,U+k) for i in (F,M) for k in range(5)]
    affected += [(U+k,j) for k in range(5) for j in [Z,M]+list(range(U,U+5))]
    for groups in range(2,len(base)):
        maxima={entry:arb(0) for entry in affected}
        for counts in multiplicities(groups):
            weight=sum(a*n for a,n in enumerate(counts,1))
            scale=rho**counts[3]*arb(input_penalty)**weight*arb(odd_penalty)**(counts[0]+counts[2])
            moments={v:min(arb(1),up(initial[v]*prod(f**n for f,n in zip(factors[v],counts)))) for v in levels}
            arbitrary=max(moments.values())
            fresh=max(up(sum((c*moments[v] for v,c in dist.items()),arb(0))/sum(dist.values()))
                      for dist in distributions.values())
            entries={(F,M):fresh/2,(M,M):arbitrary/2,(M,Z):arbitrary/(2*m)}
            for k,v in enumerate(levels):
                entries[F,U+k]=fresh*spectrum[v]/(2*m)
                entries[M,U+k]=arbitrary*spectrum[v]/(2*m)
                entries[U+k,M]=moments[v]/2
                # Uniform-class domination gives lazy cancellation <=
                # pointwise output factor / class size, independently of
                # the feedback law. Add the averaged refresh contribution.
                peak=(-arb(tilt)*abs(v-weight)).exp()
                entries[U+k,Z]=peak/(2*spectrum[v])+moments[v]/(2*m)
                for l,w in enumerate(levels):
                    entries[U+k,U+l]=moments[v]*spectrum[w]/(2*m)
            for entry,value in entries.items():
                maxima[entry]=max(maxima[entry],up(value*scale))
        for (i,j),value in maxima.items():
            base[groups][i,j]=min(base[groups][i,j],value)
    return base


if __name__=='__main__':
    self_test()
