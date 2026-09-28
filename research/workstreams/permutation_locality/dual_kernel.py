"""Exact weighted-kernel shell caps, subtracting the two known endpoint words.

Rational witnesses only. The weight factor must be nonnegative on every
remaining allowed codeword weight. No solver tolerance enters a bound.
"""
from fractions import Fraction as Q
from math import comb,log2
from flint import fmpq_mat,fmpq
from shortened_bound import krawtchouk
from joint_support import span


def rational(x):return Q(int(x.numerator),int(x.denominator))


class Kernel:
    def __init__(self,n,k,strength,allowed,factor=None):
        self.n=n;self.k=k;self.allowed=allowed
        degree=(strength-(2 if factor is not None else 0))//2
        self.degree=degree
        assert degree>=0 and all(0<w<n for w in allowed)
        self.weights=[1 if factor is None else (w-factor)*(n-factor-w) for w in range(n+1)]
        assert all(self.weights[w]>=0 for w in allowed)
        self.rows=[krawtchouk(n,w)[:degree+1] for w in range(n+1)]
        # Total codeword moment, minus the known zero and all-one words.
        gram=[]
        for i in range(degree+1):
            row=[]
            for j in range(degree+1):
                value=Q(sum(comb(n,w)*self.weights[w]*self.rows[w][i]*self.rows[w][j]
                            for w in range(n+1)),1<<(n-k))
                value-=sum(self.weights[w]*self.rows[w][i]*self.rows[w][j] for w in (0,n))
                row.append(fmpq(value.numerator,value.denominator))
            gram.append(row)
        self.gram=fmpq_mat(gram)
        self.inverse=None
        for size in range(degree+1,0,-1):
            candidate=fmpq_mat([row[:size] for row in gram[:size]])
            try: inverse=candidate.inv()
            except ZeroDivisionError: continue
            self.gram=candidate;self.inverse=inverse;self.degree=size-1
            self.rows=[row[:size] for row in self.rows]
            break

    def bound(self,w):
        if self.inverse is None or self.weights[w]<=0:return None
        vector=fmpq_mat([[x] for x in self.rows[w]])
        coefficients=self.inverse*vector
        at_w=(vector.transpose()*coefficients)[0,0]
        if at_w==0:return None
        cost=(coefficients.transpose()*self.gram*coefficients)[0,0]
        assert cost>0
        bound=rational(cost/(self.weights[w]*at_w*at_w))
        # f(x)=W(x)p(x)^2 / (W(w)p(w)^2) is nonnegative on all allowed
        # nodes. The numerator is its exact, known codeword moment.
        return bound.numerator//bound.denominator


def improve(spectrum):
    result=spectrum[:]
    allowed=list(range(30,227,2))
    for factor in (None,0,16,28,30):
        kernel=Kernel(256,128,37,allowed,factor)
        changed=0
        for w in allowed:
            bound=kernel.bound(w)
            if bound is not None and bound<result[w]:
                result[w]=bound;changed+=1
        assert result==result[::-1]
        print('Exact weighted dual kernel factor',factor,'changed',changed,'shells',flush=True)
    return result


def self_test():
    checks=0
    # Even-parity [8,7,2]: dual distance eight gives strength seven.
    # Self-dual extended Hamming [8,4,4]: strength three.
    for basis,n,strength,d in (([(1<<i)|128 for i in range(7)],8,7,2),
                              ([255,15,51,85],8,3,4)):
        code=span(basis);k=len(basis)
        for factor in (None,0,d):
            kernel=Kernel(n,k,strength,list(range(d,n-d+1,2)),factor)
            for i in range(kernel.degree+1):
                for j in range(kernel.degree+1):
                    exact=sum(kernel.weights[w.bit_count()]*kernel.rows[w.bit_count()][i]*kernel.rows[w.bit_count()][j]
                              for w in code if w not in (0,(1<<n)-1))
                    assert kernel.gram[i,j]==exact;checks+=1
            for w in kernel.allowed:
                upper=kernel.bound(w)
                if upper is not None:
                    assert sum(x.bit_count()==w for x in code)<=upper;checks+=1
    print('Weighted kernel exact code moments and shell checks:',checks,'passed',flush=True)


if __name__=='__main__':
    self_test()
    from dual_shortening import verify_bch_premise
    from dual_moments import dual_shell_caps
    verify_bch_premise()
    old=dual_shell_caps();new=improve(old)
    print('Dual shell log2 caps before/after',[(w,log2(old[w]),log2(new[w])) for w in (30,32,38,40,48,56,64,80,96,128)],flush=True)
