"""Exact integer count boxes for the unchanged layered inner.

Restricted coverage only. R and W are fixed per box; coordinates are z,r,f,e,B.
Uses monotonicity of the existing route envelope and rational local norms.
No generated reports, searches on disk, or whole-code certificate.
"""
from functools import lru_cache
from itertools import product
import time
from layered_general_components import (M,b,N,retained_failures,rational_operator,
                                        general_exact,envelopes,mixed_dimension)
from layered_component_markers import BANK,depth_of
from check_layered_component_count import log_upper_units


@lru_cache(maxsize=128)
def coefficients(R,allow_z,allow_r,bank_index,power=512):
    name,m,markers=BANK[bank_index];depth=depth_of(markers)
    size=1<<depth;contexts=size//2;P=(R*M-b)//2;KR=R*M-N//2;W=400*N//1000
    h,rows,ends=rational_operator((R,int(allow_z),int(allow_r),0,0,0,0,0),markers,m,depth)
    factors=[]
    for c in range(contexts):
        first=h[2*c][1]*markers[5]
        if allow_z or allow_r:first+=1
        factors.append((first*(1+h[2*c+1][1]*markers[5]),b//size))
    factors.extend((max(row.values()),(P-b//2)//size) for row in rows)
    factors.extend((max(row.values()),b//size) for row in ends)
    constant=power*(1+m*b+(32-m)*KR+(m-32)*P)
    constant+=sum(count*log_upper_units(value,power) for value,count in factors)
    constant+=W*log_upper_units(1/markers[5],power)
    rates=(32-m,-32,m-32,32,0)
    slopes=tuple(power*rate+log_upper_units(1/marker,power) for rate,marker in zip(rates,markers[:5]))
    return constant,slopes,m*power


def box_bound(R,lo,hi,bank_index,env=None,power=512):
    assert len(lo)==len(hi)==5 and all(0<=a<=c for a,c in zip(lo,hi))
    if env is None:env=envelopes(R)
    P=(R*M-b)//2
    # Every dimension-envelope term is decreasing in z,e and increasing in r,f.
    dim=mixed_dimension(R,P,lo[0],hi[1],hi[2],lo[3],env)
    constant,slopes,rate=coefficients(R,hi[0]>0,hi[1]>0,bank_index,power)
    value=constant+rate*dim+sum(a*(upper if a>=0 else lower) for a,lower,upper in zip(slopes,lo,hi))
    return -((-value)//power)


def ceil_log_count(lo,hi):
    count=1
    for a,c in zip(lo,hi):count*=c-a+1
    return (count-1).bit_length()


def run():
    started=time.monotonic();profiles=retained_failures();cache={};receipts=[];corner_checks=0
    for case,p in enumerate(profiles):
        R,z,r,f,e,B,_,_=p
        if R not in cache:cache[R]=envelopes(R)
        center=(z,r,f,e,B)
        lo=tuple(0 if i==0 or (i==1 and r==0) else max(0,v-128) for i,v in enumerate(center))
        hi=tuple(0 if i==0 or (i==1 and r==0) else v+128 for i,v in enumerate(center))
        best,index=min((box_bound(R,lo,hi,j,cache[R]),j) for j in range(len(BANK)))
        union=best+ceil_log_count(lo,hi)
        assert union<-100,(case,lo,hi,best)
        # The affine box majorant must dominate direct exact evaluation at
        # every corner, including boxes that keep r=0 fixed.
        _,m,markers=BANK[index]
        for point in product(*(sorted({a,c}) for a,c in zip(lo,hi))):
            zz,rr,ff,ee,bb=point;dim=mixed_dimension(R,(R*M-b)//2,zz,rr,ff,ee,cache[R])
            direct=general_exact((R,zz,rr,ff,ee,bb,dim,0),markers,m,depth_of(markers))
            assert direct<=best
            corner_checks+=1
        receipts.append((case,best,union,index))
        if time.monotonic()-started>55:raise TimeoutError('55s neighborhood/corner cap')
    # Continuous integer coverage z=0..8192, retaining positive r. The first
    # cell permits00 pairs even at its z=0 endpoint, a safe larger sum.
    p=profiles[35];R,_,r,f,e,B,_,_=p;env=cache[R];intervals=[];cursor=0
    for left in range(0,8193,256):
        right=min(8192,left+255);lo=(left,r,f,e,B);hi=(right,r,f,e,B)
        assert left==cursor
        upper=box_bound(R,lo,hi,1,env)
        assert upper<-100,(left,right,upper)
        intervals.append((left,right,upper));cursor=right+1
    assert cursor==8193
    interval_union=max(v[2] for v in intervals)+(8193-1).bit_length()
    assert interval_union<-100
    return dict(neighborhood_boxes=len(receipts),radius=128,
                worst_box=max(receipts,key=lambda x:x[2]),corner_checks=corner_checks,
                empty_count_interval=(0,8192),gapless_segments=len(intervals),
                empty_interval_log2_union_upper=interval_union,
                elapsed_seconds=time.monotonic()-started,whole_code_certificate=False)


if __name__=='__main__':print(run())
