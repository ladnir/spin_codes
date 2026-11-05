"""Common-variable switch-gradient diagnostic; FLOAT screen, not a certificate.

For positive pair polynomials F_ab(x), a matching switch is log(A/B),
where A and B each multiply two F_ab. Clearing the positive denominator
turns a squared-gradient bound into coefficientwise polynomial inequalities.
No construction change, generated files, or encoding timings.
"""
import math
import shutil
import time
from fractions import Fraction as F
from itertools import product
from random_pairing_flow_quadratic import TYPES, LAW, OUTPUT, LOSS


def switch_keys():
    keys=set()
    for a,b,x,c in product(range(5),repeat=4):
        coefficients={}
        for pair,sign in (((x,a),1),((x,b),-1),((c,a),-1),((c,b),1)):
            pair=tuple(sorted(pair));coefficients[pair]=coefficients.get(pair,0)+sign
        key=tuple(sorted((pair,k) for pair,k in coefficients.items() if k))
        if key:
            if key[0][1]<0:key=tuple((pair,-k) for pair,k in key)
            keys.add(key)
    return tuple(sorted(keys))


def multiply(left,right):
    result={}
    for a,x in left.items():
        for b,y in right.items():
            key=tuple(i+j for i,j in zip(a,b));result[key]=result.get(key,0)+x*y
    return result


def pair_polynomials(exact=False):
    result={}
    for a,b in TYPES:
        poly={}
        for j in range(3):
            weight=(LAW[j]*F(11,10)**LOSS[a,b,j]*F(6,5)**j if exact else
                    float(LAW[j])*math.exp(.1*LOSS[a,b,j]+.2*j))
            powers=OUTPUT[a,b,j];poly[powers]=poly.get(powers,0)+weight
        result[a,b]=poly
    return result


def quotient_polynomials(key,pairs):
    one={(0,)*5:1};a=one;b=one
    for pair,coefficient in key:
        for _ in range(abs(coefficient)):
            if coefficient>0:a=multiply(a,pairs[pair])
            else:b=multiply(b,pairs[pair])
    assert {sum(k) for k in a}=={sum(k) for k in b}
    return a,b


def cleared_gradient(a,b):
    # D=AB; N_i=(x_i partial_i A)B-A(x_i partial_i B).
    denominator={};numerator={}
    for pa,ca in a.items():
        for pb,cb in b.items():
            powers=tuple(i+j for i,j in zip(pa,pb));weight=ca*cb
            denominator[powers]=denominator.get(powers,0)+weight
            vector=numerator.setdefault(powers,[0]*5)
            for i in range(5):vector[i]+=weight*(pa[i]-pb[i])
    return denominator,numerator


def squared_ratio(denominator,numerator):
    # Every monomial of sum N_i^2 is also present in D^2.
    den2={};num2={};items=list(denominator.items())
    for index,(pa,ca) in enumerate(items):
        va=numerator[pa]
        for j in range(index,len(items)):
            pb,cb=items[j];vb=numerator[pb];factor=1 if j==index else 2
            powers=tuple(i+k for i,k in zip(pa,pb))
            den2[powers]=den2.get(powers,0)+factor*ca*cb
            num2[powers]=num2.get(powers,0)+factor*sum(x*y for x,y in zip(va,vb))
    assert all(x>0 for x in den2.values()) and den2.keys()==num2.keys()
    worst=max(den2,key=lambda k:num2[k]/den2[k])
    return max(0,num2[worst]/den2[worst]),worst,den2,num2


def evaluate(poly,x):
    return sum(c*math.prod(t**p for t,p in zip(x,k)) for k,c in poly.items())


def sharpness_test():
    """Unique leading monomials prove the norm-squared limit, for any positive coefficients."""
    pairs=pair_polynomials(True)
    key=(((2,2),1),((2,4),-2),((4,4),1))
    expected={(2,2):{(0,2,0,0,0),(0,0,0,2,0),(0,0,0,0,2)},
              (2,4):{(0,1,0,1,0),(0,0,0,0,2)},
              (4,4):{(0,0,0,2,0),(0,0,0,0,2)}}
    exponents=(0,3,0,-1,0);limit=[0]*5
    for pair,coefficient in key:
        assert set(pairs[pair])==expected[pair]
        order={powers:sum(a*b for a,b in zip(powers,exponents)) for powers in pairs[pair]}
        leading=[powers for powers,value in order.items() if value==min(order.values())]
        assert len(leading)==1
        for i in range(5):limit[i]+=coefficient*leading[0][i]
    assert tuple(limit)==(0,0,0,4,-4) and sum(x*x for x in limit)==32
    return dict(switch=key,positive_variable_path='x=(1,epsilon^3,1,epsilon^-1,1), epsilon -> 0+',
                limiting_log_gradient=limit,squared_norm=32,
                smaller_uniform_global_constant_ruled_out=True,code_obstruction=False)


def exact_tests():
    pairs=pair_polynomials(True);keys=switch_keys();checks=0
    # Exercise different supports and repeated factors, not just one switch.
    for key in keys[::max(1,len(keys)//8)]:
        a,b=quotient_polynomials(key,pairs);d,n=cleared_gradient(a,b)
        bound,_,d2,n2=squared_ratio(d,n)
        assert all(bound*d2[k]>=n2[k] for k in d2);checks+=1
        for x in ((F(1),)*5,tuple(F(i+1,3) for i in range(5)),tuple(F(5-i,2) for i in range(5))):
            av=evaluate(a,x);bv=evaluate(b,x);dv=evaluate(d,x)
            assert dv==av*bv
            norm=0
            for i in range(5):
                direct=evaluate({p:c*p[i] for p,c in a.items()},x)/av-evaluate({p:c*p[i] for p,c in b.items()},x)/bv
                actual=evaluate({p:v[i] for p,v in n.items()},x)/dv
                assert direct==actual;norm+=actual*actual;checks+=1
            assert norm==evaluate(n2,x)/evaluate(d2,x)<=bound;checks+=1
    return checks


def main():
    assert shutil.disk_usage('.').free>=20*(1<<30),'disk guard'
    start=time.monotonic();deadline=start+40
    checks=exact_tests();pairs=pair_polynomials();keys=switch_keys();records=[]
    for key in keys:
        if time.monotonic()>deadline:break
        a,b=quotient_polynomials(key,pairs);d,n=cleared_gradient(a,b)
        bound,powers,_,_=squared_ratio(d,n);records.append((bound,key,powers))
    print(dict(exact_identity_checks=checks,switches_checked=len(records),switches_total=len(keys),
               numerical_screen_only=True,seconds=time.monotonic()-start),flush=True)
    print(dict(generic_squared_lipschitz=32,
               proposed_max=max(x[0] for x in records),worst=sorted(records,reverse=True)[:5]),flush=True)
    print(dict(exact_leading_monomial_witness=sharpness_test()),flush=True)
    assert time.monotonic()-start<45


if __name__=='__main__':main()
