"""Exact rank checks for the systematic-inner candidate; no global certificate."""
from fractions import Fraction
from functools import lru_cache
from math import comb
import random
import unittest
from hd_sumrank_maps import construction,Q,transpose
from hd_blockwise_shared_transfer import mm,rank


@lru_cache(None)
def maps():
    m=construction(16,True,True)
    local_a=mm(m['local_f'],m['local_a'],Q)
    a=tuple(row for block in m['blocks']
            for row in mm(local_a,transpose(block),Q))
    return a,transpose(a),local_a


class IdentityFeedthroughTests(unittest.TestCase):
    def test_self_orthogonality(self):
        a,c,local_a=maps()
        self.assertEqual(mm(transpose(local_a),local_a,Q),((0,)*4,)*4)
        self.assertEqual(mm(c,a,Q),((0,)*16,)*16)
        self.assertEqual(rank(c,Q),16)
        f=tuple(tuple(int(i==j) for j in range(64)) for i in range(64))
        h=0
        for matrix in (f,a,c):
            for row in matrix:
                for x in row:h=((h^x)*0x100000001b3)&((1<<64)-1)
        self.assertEqual(h,1274355414671835662)

    def test_exact_profile_implementation(self):
        from hd_sumrank_transfer import profiles
        from hd_identity_feedthrough_screen import local_moments
        qlocal,flocal,free=local_moments(.75)
        for profile,j,mult,k,t,v in profiles(Fraction(3,4),state=4,tiles=1,feedthrough='identity'):
            self.assertEqual(mult,comb(8,j))
            self.assertEqual(t,(Fraction(3,4)+Fraction(1,4*Q))**j)
            self.assertAlmostEqual(float(k)/sum(qlocal[j,h]*Q**(-h) for h in range(5)),1,places=12)
            fv=sum(flocal[j,h]*(Q**(-h)-Q**(-4))/(1-Q**(-4)) for h in range(5))
            self.assertAlmostEqual(float(v)/fv,1,places=12)

    def test_mask_ranks(self):
        a,c,_=maps();rng=random.Random(9528)
        for _ in range(160):
            u=sorted(rng.sample(range(64),rng.randrange(65)))
            z=sorted(rng.sample(range(64),rng.randrange(65)))
            r0=len(set(u)&set(z))
            rho=min(16,sum(min(4,len({x for x in u if x//8==i}-set(z))) for i in range(8)))
            eta=min(16,sum(min(4,len({x for x in z if x//8==i}-set(u))) for i in range(8)))
            f=tuple(tuple(int(i==j) for j in u) for i in z)
            quiet=f+tuple(tuple(row[j] for j in u) for row in c)
            flat=tuple(row+a[i] for row,i in zip(f,z))
            self.assertEqual(rank(quiet,Q),r0+rho,(u,z))
            self.assertEqual(rank(flat,Q),r0+eta,(u,z))

    def test_local_moments(self):
        from hd_identity_feedthrough_screen import local_moments
        z=Fraction(3,4);quiet,flat,free=local_moments(float(z))
        _,_,a=maps();c=transpose(a)
        for j in range(9):
            u=tuple(range(j));k=Fraction(0);t=Fraction(0);v=Fraction(0)
            for mask in range(256):
                out=tuple(i for i in range(8) if mask>>i&1)
                f=tuple(tuple(int(i==x) for x in u) for i in out)
                quiet_rank=rank(f+tuple(tuple(row[x] for x in u) for row in c),Q)
                joint=rank(tuple(row+a[i] for row,i in zip(f,out)),Q)
                base=len(set(u)&set(out))
                weight=z**(8-len(out))*(1-z)**len(out)
                k+=weight*Fraction(1,Q**quiet_rank)
                t+=weight*Fraction(1,Q**base)
                v+=weight*Fraction(Q**(4-joint+base)-1,Q**base*(Q**4-1))
            qk=sum(quiet[j,h]*Q**(-h) for h in range(5))
            fv=sum(flat[j,h]*(Q**(-h)-Q**(-4))/(1-Q**(-4)) for h in range(5))
            self.assertAlmostEqual(qk/float(k),1,places=12)
            self.assertAlmostEqual(free[j]/float(t),1,places=12)
            self.assertAlmostEqual(fv/float(v),1,places=12)
            self.assertEqual(t,(z+(1-z)/Q)**j)


if __name__=='__main__':unittest.main()
