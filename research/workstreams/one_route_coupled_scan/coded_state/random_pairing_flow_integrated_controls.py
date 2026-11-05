"""Scalar degree elevation after integrating common histogram controls.

No change to the Gaussian envelope or encoder. Each simplex retains its
actual barycentric weights instead of immediately taking the largest
integrated coefficient. Standard-library, bounded, source-only research.
"""
import math
import time
from fractions import Fraction as F
from random_pairing_flow_envelope_verify import I, weak_compositions
import random_pairing_flow_envelope_verify as verifier


def elevate_once(coefficients,degree,size):
    result={}
    for beta in weak_compositions(degree+1,size):
        value=I()
        for i,count in enumerate(beta):
            if count:
                alpha=list(beta);alpha[i]-=1
                value+=coefficients[tuple(alpha)]*count
        result[beta]=value/(degree+1)
    return result


def integrated_upper(controls,vertex_count,simplices,degree,target_degree):
    """Bound a common-control polynomial on the union of its simplices."""
    assert 2<=degree<=target_degree<=12 and simplices
    assert sum(math.comb(target_degree+len(s)-1,len(s)-1) for s in simplices)<=100000
    maxima=[];count=0;original=[]
    for simplex in simplices:
        assert time.monotonic()<verifier.VERIFY_DEADLINE,'integrated control deadline'
        size=len(simplex);local={}
        for alpha in weak_compositions(degree,size):
            if degree==2:
                key=tuple(i for i,a in zip(simplex,alpha) for _ in range(a))
            else:
                key=[0]*vertex_count
                for i,a in zip(simplex,alpha):key[i]=a
                key=tuple(key)
            local[alpha]=controls[key]
        original.append(max(x.hi for x in local.values()))
        for current in range(degree,target_degree):local=elevate_once(local,current,size)
        maxima.append(max(x.hi for x in local.values()));count+=len(local)
    metadata=dict(integrated_control_degree=target_degree,integrated_scalar_controls=count)
    if max(original)>0:metadata['unelevated_simplex_log_upper']=str(I(max(original)).ln().hi)
    return I(max(maxima)),metadata


def evaluate_exact(coefficients,weights):
    degree=sum(next(iter(coefficients)));result=F(0)
    for alpha,value in coefficients.items():
        factor=F(math.factorial(degree),math.prod(math.factorial(x) for x in alpha))
        result+=value*factor*math.prod(x**a for x,a in zip(weights,alpha))
    return result


def tests():
    checks=0
    for size in (1,2,3,5):
        for degree in (2,3,4):
            exact={a:F(sum((i+1)*x*x for i,x in enumerate(a))-3,7) for a in weak_compositions(degree,size)}
            values={a:I(x) for a,x in exact.items()}
            for current in range(degree,degree+3):
                old=values;values=elevate_once(values,current,size)
                next_exact={}
                for beta in values:
                    total=F(0)
                    for i,count in enumerate(beta):
                        if count:
                            alpha=list(beta);alpha[i]-=1;total+=count*exact[tuple(alpha)]
                    next_exact[beta]=total/(current+1)
                    assert F(values[beta].lo)<=next_exact[beta]<=F(values[beta].hi);checks+=1
                for weights in ([F(1,size)]*size,[F(int(i==0)) for i in range(size)]):
                    assert evaluate_exact(exact,weights)==evaluate_exact(next_exact,weights);checks+=1
                assert max(x.hi for x in values.values())<=max(x.hi for x in old.values())+I('1e-25').hi
                exact=next_exact
    # The sole large coefficient cannot be realized with unit weight.
    # 10*2*x*(1-x) has supremum5, whereas its degree2 maximum coefficient is10.
    controls={(0,0):I(0),(0,1):I(10),(1,1):I(0)}
    upper,detail=integrated_upper(controls,2,((0,1),),2,12)
    assert upper.hi<6 and upper.lo>=5;checks+=1
    # Every simplex is included; an unused extra control need not affect it.
    controls={(i,j):I(1000 if (i,j)==(0,2) else i+j) for i in range(3) for j in range(i,3)}
    upper,_=integrated_upper(controls,3,((0,1),(1,2)),2,6)
    assert 4<=upper.hi<5;checks+=1
    return checks


if __name__=='__main__':print(dict(exact_integrated_control_checks=tests()),flush=True)
