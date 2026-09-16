"""Independent finite checks for support packing and the new Fourier cap."""
from collections import Counter
from fractions import Fraction as F
import math
import random
import unittest
from flint import arb,ctx
import audit_overlaps as audit
import overlap_transfer as transfer
from diagnose_kernel import weight_three_lower
import test_asymmetric_outward as toy
model=transfer.model


def packing(columns,s,reverse=False):
    groups=[]
    for original in reversed(columns) if reverse else columns:
        for group in groups:
            if len(group)==s:continue
            candidate=group+[original]
            if model.independent.search.rank(candidate)>len(group):
                group.append(original)
                break
        else:
            if original:groups.append([original])
    return sum(len(g)==s for g in groups)


def reference(a,b,s):
    groups=Counter()
    for dual in range(1,1<<s):
        selected=[i for i,c in enumerate(b) if (c&dual).bit_count()%2]
        complement=[i for i in range(len(a)) if i not in selected]
        left=[a[i] for i in selected];right=[a[i] for i in complement]
        d=max(packing(left,s),packing(left,s,True))
        e=max(packing(right,s),packing(right,s,True))
        for q in range(1,1<<s):
            assert sum((c&q).bit_count()%2 for c in left)>=d
            assert sum((c&q).bit_count()%2 for c in right)>=e
        groups[len(selected),d,e]+=1
    return sorted([*key,n] for key,n in groups.items())


class OverlapTests(unittest.TestCase):
    def setUp(self):ctx.prec=256

    def test_native_matches_independent_small_enumeration(self):
        rng=random.Random(813)
        for t,s in ((5,3),(8,3),(13,4),(19,5),(64,3),(128,3)):
            a=[rng.randrange(1,1<<s) for _ in range(t)]
            b=[rng.randrange(1,1<<s) for _ in range(t)]
            actual=audit.run_native(transfer.HERE/'packing_audit.exe',a,b,s)
            self.assertEqual(actual['groups'],reference(a,b,s))
            self.assertEqual(actual['states'],(1<<s)-1)

    def test_exact_toy_weighted_laws(self):
        check=toy.OutwardTests();check.setUp()
        groups=reference(check.a,check.b,3)
        for theta in (F(1,10),F(1,2),F(9,10)):
            for z in (F(1,2),F(9,10)):
                lam=-model.number(z).log()
                matrix=transfer.bernoulli(check.engine,model.number(theta),lam,groups)
                probabilities=[theta**x.bit_count()*(1-theta)**(5-x.bit_count()) for x in range(32)]
                check.check(matrix,probabilities,z)

    def test_cap_is_tighter(self):
        engine=model.Engine(20);groups=transfer.load_groups()
        for theta,lam in ((F(1,10),F(1,2)),(F(1,4),F(3,4)),(F(3,5),F(1,100))):
            old=engine.bernoulli(model.number(theta),model.number(lam))
            new=transfer.bernoulli(engine,model.number(theta),model.number(lam),groups)
            for a,b in zip(new,old):self.assertLessEqual(a,b+arb(2)**-200)

    def test_weight_three_convex_floor(self):
        engine=model.Engine(20)
        self.assertTrue(all(c.bit_count()==3 for c in engine.columns))
        for rho in (F(1,5),F(1,2),F(77,100),F(9,10)):
            # Exact rationals independently check the mean-weight interpolation.
            floor=F(0)
            for j in range(20):
                n=math.comb(19,j)
                count=128*sum(math.comb(3,h)*math.comb(16,j-h) for h in (1,3) if 0<=j-h<=16)
                degree,remainder=divmod(count,n)
                floor+=(n-remainder)*rho**degree+remainder*rho**(degree+1)
            actual=1+sum(n*rho**w for w,n in engine.b_spectrum.items())
            self.assertLessEqual(floor,actual)
            self.assertTrue(weight_three_lower(model.number(rho)).contains(model.number(floor)))


if __name__=='__main__':unittest.main()
