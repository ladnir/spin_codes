"""Exact cone checks and explicit two-step tests for the memory envelope."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import exp
import numpy as np
from flint import ctx

from group_moment import maps
from occupancy_model import local_data
from occupancy_memory import prepare,epoch_operators,self_test,Z,M,C,U,TERMINAL
from two_group_screen import log_power_moment


def cone_tests():
    # A nonzero submeasure is bounded simultaneously by mass and by density.
    # Convolution with a probability measure does not increase its density
    # cap. Removing zero preserves that statement. Nonconstant output weights
    # need only be bounded pointwise by f.
    cases=0
    n=8
    for weights in ((0,1,1,0,2,0,0,1),(0,0,1,3,0,1,0,1)):
        source=[F(x,sum(weights)) for x in weights]
        for raw in ((0,1,1,1,0,0,0,0),(1,0,2,0,0,0,1,0)):
            feedback=[F(x,sum(raw)) for x in raw]
            density=max(source)
            # Independent arbitrary output multipliers in [0,f].
            f=F(7,8)
            for seed in range(5):
                output=lambda x,y: f*F((3*x+5*y+seed)%9,8)
                lazy=[sum((source[x]*feedback[x^z]*output(x,x^z) for x in range(n)),F(0)) for z in range(n)]
                assert lazy[0]<=f*density*(1-feedback[0])
                assert sum(lazy[1:])<=f*sum(source)
                assert max(lazy[1:])<=f*density
                convolution=[sum((source[x]*feedback[x^z] for x in range(n)),F(0)) for z in range(n)]
                assert max(lazy[1:])<=f*max(convolution[1:])
                assert lazy[0]<=f*convolution[0]
                cases+=1
    # Mixtures of fresh laws need only one coefficient: max over source
    # shapes bounds cancellation and mature density without summing shapes.
    laws=[[F(int(x in supp),len(supp)) for x in range(n)] for supp in ((1,2),(3,5,7),(4,6))]
    for coefficients in product((F(0),F(1,3),F(2,3)),repeat=3):
        if not sum(coefficients):
            continue
        source=[sum(coefficients[j]*laws[j][x] for j in range(3)) for x in range(n)]
        for fresh in laws:
            conv=lambda a: [sum(a[x]*fresh[x^z] for x in range(n)) for z in range(n)]
            actual=conv(source)
            constituents=[conv(p) for p in laws]
            assert actual[0]<=sum(coefficients)*max(p[0] for p in constituents)
            assert max(actual[1:])<=sum(coefficients)*max(max(p[1:]) for p in constituents)
            cases+=1
    print('Exact density, tilted-convolution, and fresh-mixture tests passed:',cases,flush=True)


def direct_two_step(prepared):
    images,columns,_=maps()
    inputs={w:[] for w in range(1,5)}
    for j in range(32):
        for mask in range(1,16):
            q=0
            for bit in range(4):
                if mask>>bit&1:
                    q^=columns[4*j+bit]
            inputs[mask.bit_count()].append((q,mask<<(4*j)))
    cases=0
    for tilt in ('0','.0016','.01'):
        operators=epoch_operators(prepared,tilt,detailed=True)
        for a,b in product(range(1,5),repeat=2):
            size=len(inputs[a])*len(inputs[b])
            lazy=Counter()
            for (x,_),(y,word) in product(inputs[a],inputs[b]):
                lazy[x^y]+=exp(-float(tilt)*(a+(images[x]^word).bit_count()))/(2*size)
            # The true refresh branch is uniform before the feedback shift.
            # Nonzero input feedback gives exactly this zero-state mass.
            refresh=sum(lazy.values())/((1<<19)-1)
            reference=np.array([1.,0.,0.,0.,0.,0.,0.,0.,0.])
            for w in (a,b):
                t=operators[1][w,]
                reference=reference@np.array([[float(t[i,j]) for j in range(9)] for i in range(9)])
            tolerance=1e-12
            assert lazy[0]+refresh<=reference[Z]+tolerance
            assert sum(v for x,v in lazy.items() if x)<=reference[M]+tolerance
            assert max(v for x,v in lazy.items() if x)<=reference[C]+tolerance
            # Each represented uniform level gives per-state refresh density
            # at least refresh; the actual shifted law can omit some states.
            spectrum=prepared[0][0][0]
            for i,w in enumerate(sorted(spectrum)):
                assert refresh<=reference[U+i]/spectrum[w]+tolerance
            assert 2*sum(lazy.values())<=reference@TERMINAL+tolerance
            cases+=1
    print('Actual-map two-step moment, zero, mass, and density checks passed:',cases,flush=True)


def terminal_test():
    matrix=np.array([[.9,.01,.02],[.01,.8,.02],[.3,.01,.4]])
    terminal=np.array([1.,1.,0.])
    for n in (1,2,3,17,256):
        expected=(np.linalg.matrix_power(matrix,n)@terminal)[0]
        assert abs(exp(log_power_moment(matrix,n,terminal))-expected)<=1e-13
        assert log_power_moment(matrix,n)==log_power_moment(matrix,n,np.ones(3))
    print('Terminal functional and existing all-mass default checks passed',flush=True)


if __name__=='__main__':
    cone_tests()
    terminal_test()
    prepared=prepare(local_data(3))
    ctx.prec=192
    self_test(prepared)
    direct_two_step(prepared)
