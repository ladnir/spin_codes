"""Exact active-row interval coverage for ONE hard count profile.

Convexity in s, not sampled points, covers each interval. No whole-code claim.
"""
import math
import time
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from layered_outer_locations import hardest,positive_rows,rows,start,objective,M,b,N
from layered_location_markers import BANK
from check_layered_component_count import log_upper_units
from layered_two_layer_components import exact_replay


def float_curve(entry,profile):
    name,m,two,p,markers=entry;R=profile[0];x=start(markers)
    x=np.append(x,math.log2(p/(1-p)))
    rr,ee=rows(x,M,m,two);x[14:18]=rr.max(axis=1);x[18:20]=ee.max(axis=1)
    value=objective(x,profile,M,m,two)*(R*M)
    slope=(32-m)*(R-128)+R*math.log2((1-p)/p)
    return [value+slope*(s-M)+(1-R)*(math.lgamma(M+1)-math.lgamma(s+1)-math.lgamma(M-s+1))/math.log(2)
            for s in range(5901,M+1)]


@lru_cache(maxsize=32)
def exact_coeff(index,power):
    profile=hardest();R,z,r,f,e,B,dim,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e
    _,m,two,p,markers=BANK[index];h,rr,ee,zero=positive_rows(markers,m,p,two)
    factors=[]
    for c in range(2):
        factors.append(((zero+p*h[2*c][1]*markers[5])*(zero+p*h[2*c+1][1]*markers[5]),b//4))
    factors.extend((max(row),(P-b//2)//(16 if two else 4)) for row in rr)
    factors.extend((max(row),b//4) for row in ee)
    factors.extend((1/v,n) for v,n in zip(markers[:6],(z,r,f,e,B,400*N//1000)))
    factors.append((1/(1-p),R*M))
    base=power*(1+m*(V+b+dim)-32*(D+e))
    base+=sum(n*log_upper_units(v,power) for v,n in factors)
    slope=power*(32-m)*(R-128)+R*log_upper_units((1-p)/p,power)
    return base,slope


def exact_point(index,s,power=2048):
    # s>=5901: min(dbar,(R-128)s)=dbar throughout this cover.
    assert 5901<=s<=M
    base,slope=exact_coeff(index,power);R=hardest()[0]
    value=base+s*slope+(R-1)*log_upper_units(F(1,math.comb(M,s)),power)
    return -((-value)//power)


def run(verify=False):
    began=time.monotonic();profile=hardest();curves=np.array([float_curve(e,profile) for e in BANK])
    chosen=curves.argmin(axis=0);best=curves.min(axis=0);segments=[];lo=5901
    labels=[int(k) if v<-7000 else -1 for k,v in zip(chosen,best)]
    for offset in range(1,len(labels)+1):
        if offset==len(labels) or labels[offset]!=labels[offset-1]:
            hi=5900+offset;segments.append((lo,hi,labels[offset-1]));lo=hi+1
    gaps=[(lo,hi) for lo,hi,k in segments if k<0]
    result=dict(selected_profile=profile,bank_size=len(BANK),gaps=gaps,
                worst_float=float(max(best)),whole_code_certificate=False)
    if verify:
        small=exact_replay(profile)-4*64*(M-5900)+M
        verified=[]
        for lo,hi,k in segments:
            if k<0:continue
            bound=max(exact_point(k,lo),exact_point(k,hi))
            assert bound<-100,(lo,hi,k,bound)
            verified.append((lo,hi,BANK[k][0],bound))
        result.update(small_active_rows=(1,5900,small),verified_intervals=verified,
                      profile_complete=not gaps,
                      log2_good_route_class_union_upper=(max([small]+[v[3] for v in verified])+13 if not gaps else None))
    else:result['segments']=[(a,c,BANK[k][0] if k>=0 else None) for a,c,k in segments]
    result['elapsed_seconds']=time.monotonic()-began
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    print(run(args.verify))
