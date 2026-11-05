"""Exact first-layer coefficient geometry for the restricted three-state model.

Positive integer powers give the six transition polynomials, with common
denominator255^128. Decimal profiles select a candidate concavity-gap
witness only; the reported integer gap is checked by exact integer powers.
This diagnoses quota representation loss, not a code-distance failure.
Fixed profiles/bands, no search over construction parameters, no files.
Parent owns every numerical job.
"""
import os
import shutil
import threading
from math import comb, floor, log2

def _qpg_power(raw, exponent):
    value = [1]
    for _ in range(exponent):
        out = [0]*(len(value)+4)
        for j,c in enumerate(raw):
            if c:
                for i,a in enumerate(value):
                    out[i+j] += c*a
        value = out
    return value

def _qpg_gap(values, choose, logs, lower, upper):
    hull=[]
    for r in range(lower,upper+1):
        hull.append(r)
        while len(hull)>=3:
            a,b,c=hull[-3:]
            if (logs[b]-logs[a])/(b-a) >= (logs[c]-logs[b])/(c-b):
                break
            hull.pop(-2)
    best=(0,None)
    for lo,hi in zip(hull,hull[1:]):
        for r in range(lo+1,hi):
            gap=((hi-r)*logs[lo]+(r-lo)*logs[hi])/(hi-lo)-logs[r]
            if gap>best[0]: best=(gap,(lo,r,hi))
    if best[1] is None:
        return dict(band=(lower,upper),gap_display=0,exact_gap_bits=0)
    lo,r,hi=best[1];A,B,D=hi-r,r-lo,hi-lo
    g=max(0,floor(best[0])-1)
    left=pow(values[lo],A)*pow(values[hi],B)*pow(choose[r],D)
    right=pow(values[r],D)*pow(choose[lo],A)*pow(choose[hi],B)*(1<<(g*D))
    if g:
        assert left>right
    return dict(band=(lower,upper),triple=(lo,r,hi),gap_display=best[0],
                exact_gap_bits=g,exact_integer_comparison=(left>right))

def quota_profile_geometry_main():
    import resource
    if shutil.disk_usage('/tmp').free<20*(1<<30): raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(1<<30,1<<30))
    timer=threading.Timer(85,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        raw=globals().get('QE_RAW')
        if raw is None:
            from quartet_three_state_quota_envelope import QE_RAW
            raw=QE_RAW
        powers=[_qpg_power(p,128) for p in raw]
        denominator_log=128*log2(255)
        choose=[comb(512,r) for r in range(513)]
        for name,selected,hole,shift in (('O',1,0,142),('P',3,2,645),('H',5,4,645)):
            first=powers[selected]
            second=[v<<shift for v in powers[hole]]
            values=[a+b for a,b in zip(first,second)]
            assert all(v>0 for v in values)
            logs=[log2(v)-denominator_log-log2(c) for v,c in zip(values,choose)]
            dominance=[a>=b for a,b in zip(first,second)]
            changes=[r for r in range(1,513) if dominance[r]!=dominance[r-1]]
            convex=concave=0
            for r in range(1,512):
                left=values[r-1]*values[r+1]*choose[r]**2
                right=values[r]**2*choose[r-1]*choose[r+1]
                convex+=left>right;concave+=left<right
            print(dict(state=name,selected_dominance_changes=changes,
                       exact_logconvex_points=convex,exact_logconcave_points=concave,
                       boundary_log2=logs[0],
                       profiles=[(r,logs[r]) for r in (0,64,128,192,255,256,320,384,448,512)],
                       coefficients_exact=True),flush=True)
            for lo,hi in ((0,255),(256,511),(0,512)):
                print(dict(state=name,**_qpg_gap(values,choose,logs,lo,hi),
                           only_quota_profile_obstruction=True,full_code_certificate=False),flush=True)
    finally:timer.cancel()


def full_source_profile_geometry_main():
    """Full first-layer actions at O/P/H, not restricted destinations.

    phi(a,w)=2**(a/4+w) gives per-hole factors2 and32 at
    sources O and P/H respectively. Mixed hole counts are retained.
    Only three source rows are examined; subsequent actions are open.
    """
    import resource
    if shutil.disk_usage('/tmp').free<20*(1<<30): raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(1<<30,1<<30))
    timer=threading.Timer(85,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        raw=globals().get('QE_RAW')
        if raw is None:
            from quartet_three_state_quota_envelope import QE_RAW
            raw=QE_RAW
        choose=[comb(512,r) for r in range(513)]
        for name,selected,hole,factor in (('O',1,0,2),('P',3,2,32),('H',5,4,32)):
            local=[a+factor*b for a,b in zip(raw[selected],raw[hole])]
            values=_qpg_power(local,128)
            assert all(v>0 for v in values)
            logs=[log2(v)-128*log2(255)-log2(c) for v,c in zip(values,choose)]
            convex=concave=0
            for r in range(1,512):
                left=values[r-1]*values[r+1]*choose[r]**2
                right=values[r]**2*choose[r-1]*choose[r+1]
                convex+=left>right;concave+=left<right
            print(dict(full_destinations=True,source=name,local_raw=local,
                       exact_logconvex_points=convex,exact_logconcave_points=concave,
                       profiles=[(r,logs[r]) for r in (0,64,128,192,255,256,320,384,448,512)],
                       full_code_certificate=False),flush=True)
            for lo,hi in ((0,255),(256,511),(0,512)):
                print(dict(source=name,full_destinations=True,
                           **_qpg_gap(values,choose,logs,lo,hi)),flush=True)
    finally:timer.cancel()

if __name__=='__main__':
    quota_profile_geometry_main()

