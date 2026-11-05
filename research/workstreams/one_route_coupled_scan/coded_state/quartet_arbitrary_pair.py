"""Two-update moments without forgetting the intermediate age histogram.

Exact toy identities and local FLOAT bounds only. Serial Peach,180s,
<256MiB, no generated files. No construction change or distance certificate.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
import quartet_fresh_count_screen as fresh
import quartet_joint_regional as joint
from quartet_full_pair import entries,logsum
from quartet_window_transfer import kernel
from quartet_two_cohort import compositions,elementary,H,D


def multiply(p,q):
    out=defaultdict(F)
    for (i,j),a in p.items():
        for (k,l),b in q.items():out[i+k,j+l]+=a*b
    return dict(out)


def exact_moments(n,R=F(1),S=F(1),h=H,d=D):
    """All (first input count,second input count) moments R^out1 S^out2."""
    w=sum(n);b3=F(d[3],h[4]);b4=F(d[4],h[4]);core={(0,0):F(1)}
    for ni,c,e in zip(n,(1,1,1,b3,b4),(d[1]*S,d[2],d[3],d[4],d[4])):
        for _ in range(ni):core=multiply(core,{(0,0):F(1),(1,0):F(c),(1,1):c*e})
    phi=math.prod(hi**ni for hi,ni in zip(h,n));out=defaultdict(F)
    for l in range(w+1):
        second=[F(0)]*(w+1)
        for k,J,_ in entries(w,l):
            selected=[core.get((k,j),F(0)) for j in range(k+1)]
            outside=elementary((w-l,l-k),(d[0]*R,d[1]*S))
            for i,a in enumerate(selected):
                for j,b in enumerate(outside):second[i+j]+=F(J,math.comb(w,k))*a*b/phi
        for last in range(w+1):
            out[l,last]=sum(F(J,math.comb(w,k))*second[k]*(h[1]*R)**(last-k)
                            for k,J,_ in entries(w,last))
    return dict(out)


@lru_cache(maxsize=32)
def first_geometry(width,counts):
    table=np.full((len(counts),width+1),-np.inf)
    for i,l in enumerate(counts):
        for k,J,logJ in entries(width,l):table[i,k]=logJ-math.log2(math.comb(width,k))
    return table


def moment_bounds(cell,first_counts,last_counts,width=512,h=None,d=None,R=1.,S=1.):
    """Uniform over source cell; minimize each second coefficient separately."""
    h=fresh.H if h is None else np.asarray(h,dtype=float)
    d=fresh.D if d is None else np.asarray(d,dtype=float)
    b3=d[3]/h[4];b4=d[4]/h[4]
    best=np.full((len(first_counts),width+1),np.inf)
    degrees=np.arange(width+1)
    ls=np.asarray(first_counts);geom=first_geometry(width,tuple(first_counts))
    for logy in np.arange(-8.,8.001,.25):
        y=2**logy
        e=np.array([1+d[1]*S*y,1+d[2]*y,1+d[3]*y,
                    b3*(1+d[4]*y),b4*(1+d[4]*y)])
        gamma=mixed.coupled_coefficients(cell,width,d=e,h=h)
        v=d*np.array([R,S,1,1,1])*y;v=v/(1+v)**2
        terms=geom+gamma[None,:]+(width-ls[:,None])*math.log2(1+d[0]*R*y)
        terms+=(ls[:,None]-degrees[None,:])*math.log2(1+d[1]*S*y)
        peak=terms.max(axis=1)
        total=peak+np.log2(mixed.positive_exp(terms-peak[:,None],'arbitrary_pair_moment').sum(axis=1))
        variance=(width-ls)*v[0]+ls*min(v[1:])
        mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
        for theta in (.5,.75,1.,1.25):
            c=1-theta*theta/12
            mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
        best=np.minimum(best,total[:,None]-degrees[None,:]*logy+np.log2(mass[:,None]))
    result={}
    for i,l in enumerate(first_counts):
        for last in last_counts:
            result[l,last]=logsum([logJ-math.log2(math.comb(width,k))+best[i,k]
                                   +(last-k)*math.log2(h[1]*R) for k,J,logJ in entries(width,last)])
    return result


def tests():
    exact=bounded=edges=0
    for w in (4,8):
        ns=list(compositions(w,5)) if w==4 else [(2,2,2,1,1),(0,0,0,0,8),(8,0,0,0,0)]
        for n in ns:
            direct=defaultdict(F);phi=math.prod(h**v for h,v in zip(H,n))
            for (mid,l,e),p in kernel(n):
                for (out,last,f),q in kernel(mid):
                    direct[l,last,out[1],out[2]]+=p*q*(1<<(e+f))*F(math.prod(h**v for h,v in zip(H,out)),phi)
            for R,S in ((F(1),F(1)),(F(1,2),F(2)),(F(2),F(1,2))):
                actual=defaultdict(F)
                for (l,last,r,s),value in direct.items():actual[l,last]+=value*R**r*S**s
                formula=exact_moments(n,R,S)
                assert formula==dict(actual);exact+=len(formula)
                bucketed=defaultdict(F)
                for (l,last,r,s),value in direct.items():bucketed[l,last,r//2,s//2]+=value
                for (l,last,r,s),value in bucketed.items():
                    endpoint=max(R**(-j) for j in (2*r,min(last,2*r+1)))
                    endpoint*=max(S**(-j) for j in (2*s,min(last-2*r,2*s+1)))
                    assert value<=formula[l,last]*endpoint;edges+=1
                if n in ((0,0,0,0,w),(w,0,0,0,0),(2,2,2,1,1)):
                    a=w-n[0];bounds=moment_bounds((a,a,n[1],n[1],n[2],n[2]),range(w+1),range(w+1),w,H,D,float(R),float(S))
                    for key,value in actual.items():
                        assert math.log2(value)<=bounds[key]+1e-9;bounded+=1
    print(dict(exact_arbitrary_pair_moments=exact,moments_vs_uniform_float_bounds=bounded,
               exact_marked_pair_edge_envelopes=edges),flush=True)


def diagnostic():
    joint.BIN=16;joint.states.cache_clear();counts=(312,408,512)
    base,scale,labels,offsets=joint.build(counts);where={a:i for i,a in enumerate(counts)}
    for a,rl,sl in ((312,64,16),(408,192,96),(512,0,0)):
        i=where[a];boxes,lookup=joint.states(a);j=lookup[rl//16,sl//16]
        cell=(a,a,*boxes[j]);source=offsets[i]+j
        new=moment_bounds(cell,counts,counts)
        for first,last in ((312,408),(408,312),(408,408),(512,312)):
            mi=where[first];oi=where[last];ms=slice(offsets[mi],offsets[mi+1]);dest=slice(offsets[oi],offsets[oi+1])
            left=base[source,ms];right=base[ms,dest];ls=max(left);rs=np.max(right)
            total=float((left/ls)@(right/rs).sum(axis=1));assert total>0
            old=2*scale+math.log2(ls)+math.log2(rs)+math.log2(total)
            print(dict(source=cell,input_counts=[first,last],ordinary_pair_log2=old,
                       coupled_moment_log2=new[first,last],local_only=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();diagnostic();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
