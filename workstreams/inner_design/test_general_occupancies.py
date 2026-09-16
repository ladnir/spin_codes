from collections import Counter
from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np
import general_occupancies as g
import transvection as tv
import weight_memory as wm


def maps(columns,s):
    spectrum=Counter(map(int,wm.all_weights(columns,s)));del spectrum[0]
    fibers=[Counter() for _ in range(len(columns)+1)]
    for word in range(1<<len(columns)):
        syndrome=0
        for p,col in enumerate(columns):
            if word>>p&1:syndrome^=col
        fibers[word.bit_count()][syndrome]+=1
    return dict(spectrum),[row[0] for row in fibers],fibers


class GeneralTests(unittest.TestCase):
    def test_bernoulli_pointwise(self):
        columns=[1,2,3,4,5];s=3;t=len(columns);m=(1<<s)-1
        spectrum,kernel,_=maps(columns,s);levels=sorted(spectrum)
        weights=wm.all_weights(columns,s)
        for theta,z in itertools.product((F(1,10),F(1,2),F(9,10)),(F(1,2),F(9,10))):
            envelope=np.exp(g.bernoulli_epoch(spectrum,kernel,float(theta),-math.log(float(z))))
            for q in range(m+1):
                actual=[F(0) for _ in range(m+1)]
                image=sum(tv.dot(q,col)<<p for p,col in enumerate(columns))
                for x in range(1<<t):
                    syn=0
                    for p,col in enumerate(columns):
                        if x>>p&1:syn^=col
                    mass=theta**x.bit_count()*(1-theta)**(t-x.bit_count())*z**((x^image).bit_count())
                    if not q:actual[syn]+=mass;continue
                    actual[q^syn]+=mass/2
                    for y in range(1,m+1):actual[y^syn]+=mass/(2*m)
                coordinate=np.zeros(len(levels)+2);coordinate[0 if q==0 else 1]=1
                bound=coordinate@envelope
                self.assertLessEqual(float(actual[0]),bound[0]+1e-11)
                shellcaps={w:bound[i+2]/spectrum[w] for i,w in enumerate(levels)}
                residual=sum(max(0.,float(actual[y])-shellcaps[int(weights[y])]) for y in range(1,m+1))
                self.assertLessEqual(residual,bound[1]+1e-11)
                if q:
                    # A single source atom of mass 1/a_w is dominated by one
                    # shell unit; verify the shell-specific cap independently.
                    index=levels.index(int(weights[q]))+2
                    for y in range(1,m+1):
                        cap=envelope[index,levels.index(int(weights[y]))+2]/spectrum[int(weights[y])]
                        self.assertLessEqual(float(actual[y])/spectrum[int(weights[q])],cap+1e-11)

    def test_all_fiber_caps_small_maps(self):
        for s,columns in ((2,[1,2,3]),(3,[1,2,3,4,5]),(3,[1,2,3,4,5,6,7])):
            spectrum,kernel,fibers=maps(columns,s)
            caps=g.fiber_caps(len(columns),s,spectrum,kernel)
            for cap,row in zip(caps,fibers):
                self.assertLessEqual(max((v for q,v in row.items() if q),default=0),int(cap['cap']))

    def test_all_epoch_weights_pointwise(self):
        columns=[1,2,3,4,5];s=3;t=len(columns);m=(1<<s)-1
        spectrum,kernel,_=maps(columns,s);levels=sorted(spectrum)
        caps=g.fiber_caps(t,s,spectrum,kernel);low=g.low_cancellation(columns,s)
        weights=wm.all_weights(columns,s)
        initial=[[F(q==0) for q in range(m+1)]]
        coords=[[1.,0.]+[0.]*len(levels)]
        for q in range(1,m+1):
            initial.append([F(y==q) for y in range(m+1)]);coords.append([0.,1.]+[0.]*len(levels))
        for i,w in enumerate(levels):
            initial.append([F(int(weights[q])==w,spectrum[w]) for q in range(m+1)])
            coords.append([0.,0.]+[float(k==i) for k in range(len(levels))])
        for rounds,z in itertools.product((1,2),(F(1,2),F(9,10))):
            envelope=np.exp(g.epochs(spectrum,kernel,columns,caps,low,-math.log(float(z)),rounds,t))
            eps=F(1,1<<rounds)
            for j in range(t+1):
                exact=[[F(0) for _ in range(m+1)] for _ in range(m+1)]
                for support in itertools.combinations(range(t),j):
                    word=sum(1<<p for p in support);syn=0
                    for p in support:syn^=columns[p]
                    for q in range(m+1):
                        image=sum(tv.dot(q,col)<<p for p,col in enumerate(columns))
                        mass=z**((word^image).bit_count())/math.comb(t,j)
                        if not q:exact[q][syn]+=mass;continue
                        exact[q][q^syn]+=eps*mass
                        for y in range(1,m+1):exact[q][y^syn]+=(1-eps)*mass/m
                for law,coordinate in zip(initial,coords):
                    actual=[sum(law[q]*exact[q][y] for q in range(m+1)) for y in range(m+1)]
                    bound=np.array(coordinate)@envelope[j]
                    self.assertLessEqual(float(actual[0]),bound[0]+1e-11)
                    shellcaps={w:bound[i+2]/spectrum[w] for i,w in enumerate(levels)}
                    residual=sum(max(0.,float(actual[q])-shellcaps[int(weights[q])]) for q in range(1,m+1))
                    self.assertLessEqual(residual,bound[1]+1e-11)


if __name__=='__main__':unittest.main()
