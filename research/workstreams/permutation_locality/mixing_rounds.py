"""Alternative ensemble: repeated independent transvections per IMT step.

This changes the inner, not merely the proof. Existing r=1 performance
measurements do not apply to this alternative without a new implementation.
"""
from fractions import Fraction
from flint import arb
from occupancy_memory import Z,F,M,C,U,rounded
from group_rank_one_verify import up


def self_test():
    for n in (3,7):
        identity=[[Fraction(int(i==j)) for j in range(n)] for i in range(n)]
        kernel=[[(identity[i][j]+Fraction(1,n))/2 for j in range(n)] for i in range(n)]
        current=identity
        for r in range(1,5):
            current=[[sum(current[i][k]*kernel[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
            p=Fraction(1,2**r)
            assert all(current[i][j]==p*identity[i][j]+(1-p)/n for i in range(n) for j in range(n))
    print('Repeated lazy/uniform kernel identities passed exact rational checks',flush=True)


def transform(base,spectrum,tilt,rounds):
    assert rounds>=1
    if rounds==1:return base
    alpha=arb(2)**(1-rounds)  # new lazy probability / (1/2)
    beta=2-alpha           # new refresh probability / (1/2)
    result=[];levels=sorted(spectrum)
    for occupancy,old in enumerate(base):
        new=old*1
        for i in range(1,old.nrows()):
            for j in range(old.ncols()):
                new[i,j]=up(old[i,j]*(beta if U<=j<U+5 else alpha))
            if occupancy:
                if i==M or i in (9,10):
                    new[i,Z]=up(beta*old[i,Z])
                elif i==F or U<=i<U+5:
                    # Old zero bound dominates lazy + refresh. Add an upper
                    # on refresh rather than subtracting a possibly loose cap.
                    refresh=min(up(old[i,U+k]/spectrum[v]) for k,v in enumerate(levels))
                    new[i,Z]=up(alpha*old[i,Z]+(beta-alpha)*refresh)
            elif U<=i<U+5:
                v=levels[i-U]
                # Empty uniform classes retain the lazy self-loop in addition
                # to the refresh term. Its coefficient is known exactly.
                new[i,i]=up(beta*old[i,i]-(beta-alpha)*(-arb(tilt)*v).exp()/2)
        result.append(rounded(new))
    return result


if __name__=='__main__':self_test()
