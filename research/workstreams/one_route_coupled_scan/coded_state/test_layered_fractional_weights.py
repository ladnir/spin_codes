"""Exact and independent floating regressions for rational interpolation.

Floating checks do not establish certificates. Rational local operators and
directed exponent bounds establish the acceptance path.
"""
from fractions import Fraction as F
import math
import numpy as np
from layered_interpolation import LogWeight,exact_tau,exponent_bounds
from layered_outer_locations import positive_rows,start,rows,objective,hardest
from cover_layered_location_regions import BANK,POWER,location_box,region_env,M,b
from layered_general_components import mixed_dimension


def run():
    markers=BANK[1][4];p=F(99,100);equal=0
    for m in (0,1,16,24,25,31,32):
        for two in (False,True):
            assert positive_rows(markers,m,p,two)==positive_rows(markers,LogWeight(F(1<<m)),p,two)
            assert exponent_bounds(m,POWER)==(m*POWER,m*POWER)
            equal+=1
    checks=0;rowchecks=0
    for i,(_,m,two,p,markers) in enumerate(BANK):
        if not isinstance(m,LogWeight):continue
        ml,mu=exponent_bounds(m,POWER)
        assert mu%2==0 and 0<=ml<=mu<=32*POWER
        assert ml/POWER<=float(m)+1e-12<=mu/POWER+1e-12
        h,rr,ee,_=positive_rows(markers,m,p,two)
        x=np.append(start(markers),math.log2(p/(1-p)))
        fr,fe=rows(x,M,m,two)
        for actual,proposed in zip(rr,fr):
            for a,v in zip(actual,proposed):
                assert abs(math.log2(a)-v)<1e-9;rowchecks+=1
        for actual,proposed in zip(ee,fe):
            for a,v in zip(actual,proposed):
                assert abs(math.log2(a)-v)<1e-9;rowchecks+=1
        x[14:18]=fr.max(axis=1);x[18:20]=fe.max(axis=1)
        for R in (190,192,193,194):
            for delta in ((0,)*5,(128,-128,-128,128,-128)):
                point=tuple(v+d for v,d in zip(hardest()[1:6],delta))
                for refined in (False,True,'factorial'):
                    d=mixed_dimension(R,(R*M-b)//2,*point[:4],region_env(R,refined))
                    profile=(R,)+point+(d,0)
                    bound=location_box(R,point,point,i,7900,8192,refined)
                    for s in (7900,8004,8061,8132,8192):
                        direct=objective(x,profile,s,m,two)*R*M
                        assert direct<=bound+1e-5,(i,R,s,refined,direct,bound)
                        checks+=1
    for tau in (F(0),F(-1),F((1<<32)+1)):
        try:LogWeight(tau)
        except AssertionError:pass
        else:raise AssertionError('invalid interpolation scale accepted')
    print(dict(integer_weight_equivalences=equal,rational_float_rows=rowchecks,
               direct_expression_majorants=checks,invalid_scales_rejected=3))


if __name__=='__main__':run()
