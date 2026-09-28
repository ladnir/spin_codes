"""Joint support CDF tilted by the sum of the four BCH row weights.

For a>1 the weight a^-w decreases. Intersect the group-support count cap
with the total-weight CDF from the fourth power of the outer enumerator.
Differences of this upper CDF are extremal weights, not true shell counts.
"""
from itertools import product
from flint import arb,fmpz_poly


def caps(spectrum,support_counts,tilt):
    a=arb(tilt);assert a>=1
    n=len(spectrum)-1
    coefficients=list(fmpz_poly(spectrum)**4)
    coefficients += [0]*(4*n+1-len(coefficients))
    coefficients[0]-=1
    cumulative=[];total=0
    for count in coefficients:
        total+=int(count);cumulative.append(total)
    powers=[arb(1)/(a**w) for w in range(4*n+1)]
    answer=[]
    for u,limit in enumerate(support_counts):
        value=arb(0);previous=0
        for w in range(4*u+1):
            current=min(limit,cumulative[w])
            value+=(current-previous)*powers[w]
            previous=current
        answer.append(int(value.upper().ceil().unique_fmpz()))
    assert all(x<=y for x,y in zip(answer,answer[1:]))
    return answer


def self_test():
    for basis in ((1,2),(3,5),(7,25,42)):
        words=[0]
        for b in basis:words+=[x^b for x in words]
        n=max(words).bit_length()
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        tuples=[xs for xs in product(words,repeat=4) if any(xs)]
        data=[((xs[0]|xs[1]|xs[2]|xs[3]).bit_count(),sum(x.bit_count() for x in xs)) for xs in tuples]
        counts=[sum(u<=i for u,w in data) for i in range(n+1)]
        for tilt in ('1','1.02','1.1','2'):
            upper=caps(spectrum,counts,tilt)
            for u in range(n+1):
                actual=sum((arb(tilt)**(-w) for v,w in data if v<=u),arb(0))
                assert actual<=upper[u] or (actual-upper[u]).contains(0)
    print('Total-weight tilted support CDF: exhaustive small-code checks passed',flush=True)


if __name__=='__main__':self_test()
