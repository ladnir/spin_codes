from fractions import Fraction as Q
from math import comb,factorial,log
from itertools import product
import unittest

import shuffle_comparison as sc


def profiles(n,k):
    if k==1:
        yield (n,)
    else:
        for a in range(n+1):
            for tail in profiles(n-a,k-1):yield (a,*tail)


def multinomial_mass(counts,law):
    coefficient=factorial(sum(counts));mass=Q(1)
    for n,p in zip(counts,law):coefficient//=factorial(n);mass*=p**n
    return coefficient*mass


class ShuffleComparisonTests(unittest.TestCase):
    def test_binomial_maximum_exhaustively(self):
        for slots in range(1,13):
            for q in range(slots+1):
                for p in (Q(0),Q(1,10),Q(1,3),Q(1,2),Q(9,10),Q(1)):
                    ratios=[];reference=Q(q,slots)*p
                    for a in range(q+1):
                        numerator=comb(q,a)*p**a*(1-p)**(q-a)
                        if not numerator:continue
                        denominator=comb(slots,a)*reference**a*(1-reference)**(slots-a)
                        ratios.append(numerator/denominator)
                    bound=sc.binomial_loss(slots,q,p)
                    self.assertEqual(max(ratios),bound)
                    self.assertAlmostEqual(sc.binomial_logloss(slots,q,float(p)),log(float(bound)),places=11)

    def test_full_category_profiles_including_fixed_packet_types(self):
        tilt=(Q(1),Q(1,3),Q(1,5),Q(1,7),Q(2))
        examples=[((Q(3,5),Q(2,5),Q(0),Q(0),Q(0)),()),
                  (tuple(Q(comb(4,j),16) for j in range(5)),()),
                  ((Q(1,3),Q(1,3),Q(1,3),Q(0),Q(0)),((4,1),)),
                  ((Q(1,3),Q(1,3),Q(1,3),Q(0),Q(0)),((3,1),(4,1)))]
        for law,fixed in examples:
            slots=7;q=4;normalizer=sum(p*t for p,t in zip(law,tilt))
            tilted=tuple(p*t/normalizer for p,t in zip(law,tilt))
            zero_slots=slots-q-sum(n for _,n in fixed)
            mean=[Q(q,slots)*p for p in tilted];mean[0]+=Q(zero_slots,slots)
            for j,n in fixed:mean[j]+=Q(n,slots)
            ratios=[]
            for active in profiles(q,5):
                numerator=multinomial_mass(active,tilted)
                if not numerator:continue
                total=list(active);total[0]+=zero_slots
                for j,n in fixed:total[j]+=n
                ratios.append(numerator/multinomial_mass(total,mean))
            geometry=(slots,q,law,fixed);bound=sc.exact(geometry,tilt)
            self.assertEqual(max(ratios),bound)
            self.assertAlmostEqual(sc.floating(geometry,list(map(float,tilt))),log(float(bound)),places=11)

    def test_incompatible_compositions_are_rejected(self):
        zero=('zero',Q(1),(Q(1),0,0,0,0),0)
        a=('a',Q(1),(Q(1,2),Q(1,2),0,0,0),1)
        b=('b',Q(1),(Q(1,3),Q(2,3),0,0,0),1)
        constant=('c',Q(1),(0,Q(1),0,0,0),1)
        with self.assertRaises(ValueError):sc.structure([zero,a,b],[5,1,1])
        with self.assertRaises(ValueError):sc.structure([zero,a,constant],[5,1,1])

    def test_capped_multinomial_mode_exhaustively(self):
        for law in ((Q(1,2),Q(1,3),Q(1,6)),(Q(0),Q(2,5),Q(3,5)),(Q(1,3),)*3):
            for caps in product(range(4),repeat=3):
                for total in range(sum(caps)+1):
                    values=[multinomial_mass(a,law) for a in profiles(total,3)
                            if all(x<=c for x,c in zip(a,caps))]
                    mode=sc.multinomial_mode(total,law,caps)
                    mass=Q(0) if mode is None else multinomial_mass(mode,law)
                    self.assertEqual(mass,max(values))

    def test_posterior_bound_for_heterogeneous_laws(self):
        laws=((Q(1),Q(0),Q(0),Q(0),Q(0)),
              (Q(3,5),Q(2,5),Q(0),Q(0),Q(0)),
              (Q(1,3),Q(1,3),Q(1,3),Q(0),Q(0)),
              tuple(Q(comb(4,j),16) for j in range(5)),
              (Q(0),Q(0),Q(0),Q(0),Q(1)))
        components=[(str(i),Q(1),law,int(i>0)) for i,law in enumerate(laws)]
        for counts in ((3,2,2,0,0),(3,1,1,1,1),(0,2,2,1,1),(0,0,0,0,7),(7,0,0,0,0)):
            geometry=sc.posterior_structure(components,list(counts))
            for tilt in ((Q(1),)*5,(Q(1),Q(1,3),Q(1,5),Q(1,7),Q(2))):
                tilted=[tuple(p*t/sum(a*b for a,b in zip(law,tilt)) for p,t in zip(law,tilt)) for law in laws]
                distribution={(0,)*5:Q(1)}
                for n,law in zip(counts,tilted):
                    for _ in range(n):
                        nxt={}
                        for a,mass in distribution.items():
                            for j,p in enumerate(law):
                                if not p:continue
                                b=list(a);b[j]+=1;b=tuple(b)
                                nxt[b]=nxt.get(b,Q(0))+mass*p
                        distribution=nxt
                mean=[sum(n*law[j] for n,law in zip(counts,tilted))/sum(counts) for j in range(5)]
                maximum=max(mass/multinomial_mass(a,mean) for a,mass in distribution.items())
                bound=sc.posterior_exact(geometry,tilt)
                self.assertGreaterEqual(bound,maximum)
                self.assertAlmostEqual(sc.posterior_floating(geometry,list(map(float,tilt))),log(float(bound)),places=11)
        # Equal packet laws must not pay for artificial hidden labels.
        a=components[1]
        one=sc.posterior_structure([a],[7])
        split=sc.posterior_structure([a,a],[3,4])
        self.assertEqual(one,split)


if __name__=='__main__':unittest.main()
