"""Exact small joint-cohort identities and a bounded full-width local screen.

Peach only, serial90s,<64MiB,no files. No end-to-end distance certificate.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
from quartet_window_transfer import kernel,gate_counts

H=(1,1,4,16,32);D=(1,4,16,64,128)


def compositions(total,parts):
    if parts==1:
        yield (total,);return
    for first in range(total+1):
        for rest in compositions(total-first,parts-1):yield (first,)+rest


def elementary(counts,ds):
    coeff=[1]
    for n,d in zip(counts,ds):
        for _ in range(n):
            after=[0]*(len(coeff)+1)
            for j,c in enumerate(coeff):after[j]+=c;after[j+1]+=d*c
            coeff=after
    return coeff


def joint_exact(counts,h=H,ds=D):
    width=sum(counts);a=width-counts[0];r=counts[1]
    old=elementary(counts[2:],ds[2:]);phi=math.prod(hi**n for hi,n in zip(h,counts))
    result=defaultdict(F)
    for (l,k),multiplicity in gate_counts(width).items():
        for v0 in range(max(0,k-a),min(k,width-a)+1):
            for v1 in range(max(0,k-v0-(a-r)),min(r,k-v0)+1):
                u=k-v0-v1
                numerator=multiplicity*math.comb(width-a,v0)*math.comb(r,v1)
                numerator*=ds[0]**v0*ds[1]**v1*h[1]**(l-k)*old[u]
                result[l,l-k+v0,v1]+=F(numerator,math.comb(width,k)*phi)
    return dict(result)


def factored_joint_exact(counts):
    width=sum(counts);a=width-counts[0];r=counts[1];allowed=defaultdict(F)
    for (l,k),multiplicity in gate_counts(width).items():
        for v0 in range(max(0,k-a),min(k,width-a)+1):
            allowed[l,k-v0]+=F(multiplicity*math.comb(width-a,v0)*D[0]**v0*H[1]**(l-k),math.comb(width,k))
    old=elementary(counts[2:],D[2:]);phi=math.prod(h**n for h,n in zip(H,counts));result={}
    for (l,u),A in allowed.items():
        for v1 in range(max(0,u-(a-r)),min(r,u)+1):
            result[l,l-u,v1]=A*F(math.comb(r,v1)*D[1]**v1*old[u-v1],phi)
    return result


@lru_cache(maxsize=640)
def joint_upper(width,a,r,s):
    old=a-r-s;result=defaultdict(F)
    for (l,k),multiplicity in gate_counts(width).items():
        for v0 in range(max(0,k-a),min(k,width-a)+1):
            for v1 in range(max(0,k-v0-(a-r)),min(r,k-v0)+1):
                degree=k-v0-v1;candidates=[]
                for z in (F(1,128),F(1,32),F(1,8),F(1,2),F(2)):
                    second=F(1+D[2]*z,H[2])**s
                    older=max(F(1+D[i]*z,H[i]) for i in (3,4))**old
                    candidates.append(second*older/z**degree)
                numerator=multiplicity*math.comb(width-a,v0)*math.comb(r,v1)
                numerator*=D[0]**v0*D[1]**v1*H[1]**(l-k)
                factor=F(numerator,math.comb(width,k)*H[1]**r)
                result[l,l-k+v0,v1]+=factor*min(candidates)
    return dict(result)


def tests():
    exact_checked=upper_checked=normalizations=0
    for width in (4,8):
        for counts in compositions(width,5):
            direct=defaultdict(F);phi=math.prod(h**n for h,n in zip(H,counts))
            for (out,l,exponent),weight in kernel(counts):
                direct[l,out[1],out[2]]+=weight*(1<<exponent)*F(math.prod(h**n for h,n in zip(H,out)),phi)
            actual=joint_exact(counts)
            assert actual==dict(direct);exact_checked+=len(actual)
            assert actual==factored_joint_exact(counts)
            a=width-counts[0];r=counts[1];s=counts[2]
            upper=joint_upper(width,a,r,s)
            assert all(value<=upper[key] for key,value in actual.items());upper_checked+=len(actual)
            unweighted=joint_exact(counts,(1,)*5,(1,)*5)
            assert sum(unweighted.values())==1<<width;normalizations+=1
    print(dict(exact_joint_cohort_coefficients=exact_checked,joint_cauchy_upper_checks=upper_checked,
               unweighted_histogram_normalizations=normalizations,all_width4_and8_histograms=True,
               whole_code_distance=False),flush=True)
    polygon_tests()


def polygon_tests():
    checked=empty=0
    for a in range(9):
        for r in range(a+1):
            for l in range(9):
                for xlo in range(0,l+1,2):
                    xhi=min(l,xlo+1)
                    for slo in range(0,l+1,2):
                        shi=min(l,slo+1)
                        actual=[(x,s) for x in range(xlo,xhi+1) for s in range(slo,shi+1)
                                if 0<=s<=r and max(0,l-x-a+r)<=s<=l-x]
                        vertices=rectangle_vertices(a,r,l,xlo,xhi,slo,shi)
                        assert bool(vertices)==bool(actual)
                        if not actual:empty+=1;continue
                        for c,d in ((1,1),(1,-1),(-1,2),(-2,-1)):
                            assert max(c*x+d*s for x,s in actual)==max(c*x+d*s for x,s in vertices)
                        checked+=1
    print(dict(exact_joint_polygon_boxes=checked,empty_boxes_checked=empty,
               factored_joint_coefficients_equal_direct=True),flush=True)


def local_envelope(a,r,s,l):
    """FLOAT one-update bound; s=None discards the second age count."""
    import numpy as np
    from quartet_history_screen import occupation,column_perron,logchoose
    logb=4;lam=707/2560;tau=381/512;h=np.array(column_perron(logb,lam,tau))
    ds=np.array((h[1],h[2],h[3],h[4]*2**(logb*(3*lam-tau)),h[4]*2**(logb*lam)))
    entries=occupation(l);ks=np.array([k for k,J,logJ in entries],dtype=float)
    known=np.array((512-a,r,0 if s is None else s))
    permitted=np.array((2,3,4) if s is None else (3,4));older=a-r-(0 if s is None else s)
    lo=np.full(ks.shape,-80.);hi=np.full(ks.shape,80.)
    for _ in range(55):
        x=(lo+hi)/2;z=np.exp2(x)
        logs=np.log2((1+ds[permitted,None]*z)/h[permitted,None])
        selected=permitted[np.argmax(logs,axis=0)]
        derivative=sum(n*d*z/(1+d*z) for n,d in zip(known,ds[:3]))
        derivative+=older*ds[selected]*z/(1+ds[selected]*z)-ks
        lo=np.where(derivative<0,x,lo);hi=np.where(derivative>=0,x,hi)
    x=(lo+hi)/2;z=np.exp2(x)
    value=sum(n*np.log2((1+d*z)/hh) for n,d,hh in zip(known,ds[:3],h[:3]))
    value+=older*np.max(np.log2((1+ds[permitted,None]*z)/h[permitted,None]),axis=0)-ks*x
    variance=sum(n*d*z/(1+d*z)**2 for n,d in zip(known,ds[:3]))
    variance+=older*np.min(ds[permitted,None]*z/(1+ds[permitted,None]*z)**2,axis=0)
    mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
    for theta in (.5,.75,1.,1.25):
        c=1-theta*theta/12
        mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    logs=value+np.log2(mass)+np.array([logJ-logchoose(512,k)+(l-k)*math.log2(h[1]) for k,J,logJ in entries])
    peak=float(np.max(logs))
    return peak+math.log2(float(np.sum(np.exp2(logs-peak))))


def local_screen():
    for a in (392,408,424):
        for r in (40,160,240):
            old=local_envelope(a,r,None,a)
            second_counts=[s for s in (0,32,64,96) if s<=a-r]
            gains={s:old-local_envelope(a,r,s,a) for s in second_counts}
            print(dict(input_count=a,fresh_count=r,old_log2_local_upper=old,
                       saved_log2_by_fixed_second_age=gains,
                       incoming_histogram_constraint_only=True,
                       propagation_not_yet_bounded=True,whole_code_distance=False),flush=True)


def cohort_vertices(a,r,l,left,right,constrained=True):
    """After previous (a,r), retain the forced range of the new type2 count."""
    if constrained:left=max(left,l-a)
    if left>right:return []
    points={left,right}
    if constrained:
        points.update(x for x in (l-r,l-a+r) if left<=x<=right)
    vertices=set()
    for x in points:
        lo=max(0,l-x-a+r) if constrained else 0
        hi=min(r,l-x) if constrained else l-x
        if lo<=hi:vertices.update(((x,lo),(x,hi)))
    return sorted(vertices)


def vertex_envelope(l,next_l,vertices):
    """Cauchy bound uniform on the joint (n1,n2) polygon, FLOAT only."""
    import numpy as np
    import quartet_fresh_count_screen as fresh
    from quartet_history_screen import occupation,logchoose
    h=fresh.H;ds=fresh.D
    pairs=np.array(vertices,dtype=float);r=pairs[:,0];s=pairs[:,1];older=l-r-s
    entries=occupation(next_l);ks=np.array([k for k,J,logJ in entries],dtype=float)
    lo=np.full(ks.shape,-80.);hi=np.full(ks.shape,80.)
    for _ in range(55):
        x=(lo+hi)/2;z=np.exp2(x)
        logs=np.log2((1+ds[:,None]*z)/h[:,None])
        selected=np.argmax(logs[3:],axis=0)+3
        old=np.max(logs[3:],axis=0)
        values=r[:,None]*logs[1]+s[:,None]*logs[2]+older[:,None]*old
        j=np.argmax(values,axis=0)
        ps=ds[:,None]*z/(1+ds[:,None]*z)
        derivative=(512-l)*ps[0]+r[j]*ps[1]+s[j]*ps[2]+older[j]*ps[selected,np.arange(len(ks))]-ks
        lo=np.where(derivative<0,x,lo);hi=np.where(derivative>=0,x,hi)
    x=(lo+hi)/2;z=np.exp2(x)
    logs=np.log2((1+ds[:,None]*z)/h[:,None]);old=np.max(logs[3:],axis=0)
    value=(512-l)*logs[0]+np.max(r[:,None]*logs[1]+s[:,None]*logs[2]+older[:,None]*old,axis=0)-ks*x
    vs=ds[:,None]*z/(1+ds[:,None]*z)**2
    variance=(512-l)*vs[0]+np.min(r[:,None]*vs[1]+s[:,None]*vs[2]+older[:,None]*np.min(vs[3:],axis=0),axis=0)
    mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
    for theta in (.5,.75,1.,1.25):
        c=1-theta*theta/12
        mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    logs=value+np.log2(mass)+np.array([logJ-logchoose(512,k)+(next_l-k)*math.log2(h[1]) for k,J,logJ in entries])
    peak=float(np.max(logs));return peak+math.log2(float(np.sum(np.exp2(logs-peak))))


def first_row(a,r,l):
    import numpy as np
    import quartet_fresh_count_screen as fresh
    from quartet_history_screen import occupation,logchoose
    coeff=fresh.active_coefficient(a,r,r,np.arange(a+1,dtype=float))
    row=np.full(len(fresh.buckets(l)),-np.inf)
    for k,J,logJ in occupation(l):
        for v0 in range(max(0,k-a),min(k,512-a)+1):
            base=logJ-logchoose(512,k)+logchoose(512-a,v0)
            base+=v0*math.log2(fresh.D[0])+(l-k)*math.log2(fresh.H[1])+coeff[k-v0]
            j=(l-k+v0)//fresh.BIN;row[j]=np.logaddexp2(row[j],base)
    return row


def two_update_screen():
    import numpy as np
    import quartet_fresh_count_screen as fresh
    fresh.BIN=8
    for l,next_l in ((392,392),(392,424),(424,392),(424,424)):
        for r in (40,160,240):
            row=first_row(408,r,l);bounds=[]
            for constrained in (False,True):
                logs=row.copy()
                for j,(left,right) in enumerate(fresh.buckets(l)):
                    vertices=cohort_vertices(408,r,l,left,right,constrained)
                    if not vertices:
                        assert not np.isfinite(row[j]), 'Cannot discard a positive incoming contribution'
                    elif np.isfinite(row[j]):logs[j]+=vertex_envelope(l,next_l,vertices)
                peak=float(np.max(logs));bounds.append(peak+math.log2(float(np.sum(np.exp2(logs-peak)))))
            print(dict(input_counts=[408,l,next_l],initial_fresh_count=r,
                       float_log2_two_update_upper=bounds[1],causal_cohort_saving=bounds[0]-bounds[1],
                       intermediate_fresh_count_summed=True,second_age_chronology_enforced=True,
                       directed_certificate=False,whole_code_distance=False),flush=True)


def rectangle_vertices(a,r,l,xlo,xhi,slo,shi):
    """Rectangle intersected with l-a+r <= n1+n2 <= l and n2<=r."""
    xlo=max(xlo,l-a);shi=min(shi,r);lower=l-a+r
    if xlo>xhi or slo>shi:return []
    points={(x,s) for x in (xlo,xhi) for s in (slo,shi)}
    for boundary in (lower,l):
        points.update((x,boundary-x) for x in (xlo,xhi))
        points.update((boundary-s,s) for s in (slo,shi))
    return sorted((x,s) for x,s in points if xlo<=x<=xhi and slo<=s<=shi and lower<=x+s<=l)


def older_coefficient(count):
    """Uniform active-old coefficient, with no incoming age2 assumption."""
    import numpy as np
    import quartet_fresh_count_screen as fresh
    assert count>=0
    if count==0:return np.zeros(1)
    ds=fresh.D[2:];h=fresh.H[2:];ks=np.arange(count+1,dtype=float)
    lo=np.full(ks.shape,-80.);hi=np.full(ks.shape,80.)
    for _ in range(55):
        x=(lo+hi)/2;z=np.exp2(x);logs=np.log2((1+ds[:,None]*z)/h[:,None])
        selected=np.argmax(logs,axis=0);d=ds[selected]
        derivative=count*d*z/(1+d*z)-ks
        lo=np.where(derivative<0,x,lo);hi=np.where(derivative>=0,x,hi)
    x=(lo+hi)/2;z=np.exp2(x)
    value=count*np.max(np.log2((1+ds[:,None]*z)/h[:,None]),axis=0)-ks*x
    variance=count*np.min(ds[:,None]*z/(1+ds[:,None]*z)**2,axis=0)
    mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
    for theta in (.5,.75,1.,1.25):
        c=1-theta*theta/12
        mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    return value+np.log2(mass)


def first_joint_row(a,r,l):
    """Sum the selected-fresh count, rather than keeping only its range."""
    import numpy as np
    import quartet_fresh_count_screen as fresh
    from quartet_history_screen import occupation,logchoose
    old=older_coefficient(a-r);allowed=np.full(min(a,l)+1,-np.inf)
    for k,J,logJ in occupation(l):
        for v0 in range(max(0,k-a),min(k,512-a)+1):
            u=k-v0
            base=logJ-logchoose(512,k)+logchoose(512-a,v0)
            base+=v0*math.log2(fresh.D[0])+(l-k)*math.log2(fresh.H[1])
            allowed[u]=np.logaddexp2(allowed[u],base)
    result=np.full((len(fresh.buckets(l)),r//fresh.BIN+1),-np.inf)
    for u,A in enumerate(allowed):
        if not math.isfinite(A):continue
        for v1 in range(max(0,u-(a-r)),min(r,u)+1):
            value=A+logchoose(r,v1)+v1*math.log2(fresh.D[1])-r*math.log2(fresh.H[1])+old[u-v1]
            i=(l-u)//fresh.BIN;j=v1//fresh.BIN
            result[i,j]=np.logaddexp2(result[i,j],value)
    return result


def joint_probability_screen():
    import numpy as np
    import quartet_fresh_count_screen as fresh
    fresh.BIN=8
    for r,l,next_l in ((40,424,424),(160,424,424),(160,424,392)):
        row=first_joint_row(408,r,l);weighted=[];control=[];cells=0
        for i,(left,right) in enumerate(fresh.buckets(l)):
            relaxed=vertex_envelope(l,next_l,cohort_vertices(408,r,l,left,right,False))
            for j in range(row.shape[1]):
                if not np.isfinite(row[i,j]):continue
                vertices=rectangle_vertices(408,r,l,left,right,8*j,8*j+7)
                assert vertices,'Positive joint contribution cannot have an empty constraint set'
                weighted.append(row[i,j]+vertex_envelope(l,next_l,vertices))
                control.append(row[i,j]+relaxed);cells+=1
        def total(values):
            peak=max(values);return float(peak+math.log2(sum(2**(x-peak) for x in values)))
        refined=total(weighted);baseline=total(control)
        print(dict(input_counts=[408,l,next_l],initial_fresh_count=r,
                   joint_cohort_cells=cells,float_log2_joint_two_update_upper=refined,
                   joint_cohort_saving_against_same_first_step=baseline-refined,
                   actual_selected_fresh_count_summed=True,
                   directed_certificate=False,whole_code_distance=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        if '--local-screen' in sys.argv:local_screen()
        if '--two-update' in sys.argv:two_update_screen()
        if '--joint-probability' in sys.argv:joint_probability_screen()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
