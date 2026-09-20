from collections import Counter
from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np
import search_feedback as search
import screen_dense
g=search.g


class AsymmetricTests(unittest.TestCase):
    def test_independent_bernoulli_transfer(self):
        a=[1,2,3,4,5];b=[3,5,6,7,1];s=3;t=5;m=7
        weights=search.wm.all_weights(a,s)
        spectrum=Counter(map(int,weights));del spectrum[0]
        dual=Counter(map(int,search.wm.all_weights(b,s)));del dual[0]
        kernel=[0]*(t+1)
        for x in range(1<<t):kernel[x.bit_count()]+=search.inject(b,x)==0
        levels=sorted(spectrum)
        for theta,z in itertools.product((F(1,10),F(1,2),F(9,10)),(F(1,2),F(9,10))):
            envelope=np.exp(screen_dense.bernoulli(spectrum,dual,kernel,float(theta),-math.log(float(z))))
            exact=[[F(0) for _ in range(m+1)] for _ in range(m+1)]
            for x in range(1<<t):
                syn=search.inject(b,x);p=theta**x.bit_count()*(1-theta)**(t-x.bit_count())
                for q in range(m+1):
                    mass=p*z**((x^search.image(a,q)).bit_count())
                    if not q:exact[q][syn]+=mass;continue
                    exact[q][q^syn]+=mass/2
                    for y in range(1,m+1):exact[q][y^syn]+=mass/(2*m)
            # Arbitrary source atoms test D; the maximal pointwise shell law
            # tests each C_v, including simultaneous mass at every source.
            laws=[([F(y==q) for y in range(m+1)],0 if q==0 else 1) for q in range(m+1)]
            laws += [([F(int(weights[q])==v,spectrum[v]) for q in range(m+1)],i+2) for i,v in enumerate(levels)]
            for law,index in laws:
                actual=[sum(law[q]*exact[q][y] for q in range(m+1)) for y in range(m+1)]
                bound=envelope[index]
                self.assertLessEqual(float(actual[0]),bound[0]+1e-11)
                caps={w:bound[i+2]/spectrum[w] for i,w in enumerate(levels)}
                residual=sum(max(0.,float(actual[q])-caps[int(weights[q])]) for q in range(1,m+1))
                self.assertLessEqual(residual,bound[1]+1e-11)

    def test_pair_counts(self):
        for columns in ([1,2,3,4,5],[1,3,5,7,9,11],[3,5,6,9,10,12]):
            counts=search.low_kernel(columns)
            for j in (3,4):
                exact=0
                for support in itertools.combinations(range(len(columns)),j):
                    x=sum(1<<p for p in support)
                    exact+=search.inject(columns,x)==0
                self.assertEqual(exact,counts[f'weight{j}'])

    def test_adjoint(self):
        a=[1,2,3,4,5];b=[3,5,6,7,1];t=5;s=3
        mixers=[(1,2),(3,3),(5,5)]
        self.assertTrue(any(search.inject(b,search.image(a,1<<j)) for j in range(s)))
        def run(words):
            q=0;out=[]
            for x,(u,v) in zip(words,mixers):
                out.append(x^search.image(a,q))
                q=g.tv.update(q,u,v)^search.inject(b,x)
            return out
        def transpose(words,wrong=False):
            r=0;out=[]
            for x,(u,v) in zip(words[::-1],mixers[::-1]):
                y=x^search.image(b,r);out.append(y)
                r=g.tv.update(r,v,u)^search.inject(a,y if wrong else x)
            return out[::-1]
        differences=0
        for i,j in itertools.product(range(3*t),repeat=2):
            x=[0]*3;u=[0]*3;x[i//t]=1<<(i%t);u[j//t]=1<<(j%t)
            lhs=sum((v&w).bit_count() for v,w in zip(run(x),u))&1
            rhs=sum((v&w).bit_count() for v,w in zip(x,transpose(u)))&1
            self.assertEqual(lhs,rhs)
            differences+=transpose(u)!=transpose(u,True)
        self.assertGreater(differences,0,'symmetric shortcut should fail for this asymmetric example')

    def test_general_weighted_transfer(self):
        a=[1,2,3,4,5];b=[3,5,6,7,1];s=3;t=5;m=7
        weights=search.wm.all_weights(a,s)
        spectrum=Counter(map(int,weights));del spectrum[0]
        dual=Counter(map(int,search.wm.all_weights(b,s)));del dual[0]
        kernel=[0]*(t+1)
        for x in range(1<<t):kernel[x.bit_count()]+=search.inject(b,x)==0
        caps=g.fiber_caps(t,s,dual,kernel);low=search.low_cancellation(a,b,s)
        levels=sorted(spectrum)
        initial=[[F(q==0) for q in range(m+1)]]
        coordinates=[[1.,0.]+[0.]*len(levels)]
        for q in range(1,m+1):
            initial.append([F(y==q) for y in range(m+1)])
            coordinates.append([0.,1.]+[0.]*len(levels))
        for i,w in enumerate(levels):
            initial.append([F(int(weights[q])==w,spectrum[w]) for q in range(m+1)])
            coordinates.append([0.,0.]+[float(i==j) for j in range(len(levels))])
        for z in (F(1,2),F(9,10)):
            envelope=np.exp(g.epochs(spectrum,kernel,b,caps,low,-math.log(float(z)),1,t))
            for j in range(t+1):
                exact=[[F(0) for _ in range(m+1)] for _ in range(m+1)]
                for support in itertools.combinations(range(t),j):
                    x=sum(1<<p for p in support);syn=search.inject(b,x)
                    for q in range(m+1):
                        mass=z**((x^search.image(a,q)).bit_count())/math.comb(t,j)
                        if not q:exact[q][syn]+=mass;continue
                        exact[q][q^syn]+=mass/2
                        for y in range(1,m+1):exact[q][y^syn]+=mass/(2*m)
                for law,coordinate in zip(initial,coordinates):
                    actual=[sum(law[q]*exact[q][y] for q in range(m+1)) for y in range(m+1)]
                    bound=np.array(coordinate)@envelope[j]
                    self.assertLessEqual(float(actual[0]),bound[0]+1e-11)
                    shellcaps={w:bound[i+2]/spectrum[w] for i,w in enumerate(levels)}
                    residual=sum(max(0.,float(actual[q])-shellcaps[int(weights[q])]) for q in range(1,m+1))
                    self.assertLessEqual(residual,bound[1]+1e-11)


if __name__=='__main__':unittest.main()
