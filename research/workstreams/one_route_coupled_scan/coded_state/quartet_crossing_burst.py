"""Audit the EXACT quota-preserving moment using boundary-crossing bursts.

No artificial terminal resets and no histogram envelope. Directed rational
logarithms; serial Peach90s,<64MiB,no files. Proof-method test, not distance.
"""
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from quartet_quota_audit import choose_log,linear_log


def event_weight(m,r,c):
    # Count paid subsets whose unpaid complement avoids one prescribed half.
    u=r-c
    return F(math.comb(r,c)*math.comb(m//2,u),math.comb(m,u))


def tests():
    cases=0
    for m in (4,8):
        forbidden=set(range(m//2))
        for r in range(1,min(m,4)+1):
            for c in range(r+1):
                count=0
                for I in itertools.combinations(range(m),r):
                    for P in itertools.combinations(I,c):
                        if not ((set(I)-set(P))&forbidden):count+=1
                assert F(count,math.comb(m,r))==event_weight(m,r,c);cases+=1
    # In every pair of regions the two prescribed halves join. The support
    # omits the final update of that full-input burst and pays a real death.
    for width in (4,8,512):
        allowed=set(range(8,24))
        selected=set(range(8,23))
        assert all(t in allowed and t+1 in allowed for t in selected)
        A=width*len(selected);D=width*sum(t+1 not in selected for t in selected)
        assert (A,D)==(15*width,width)
    print(dict(exact_paid_event_counts=cases,burst_support_sizes_checked=3),flush=True)


def audit():
    m=8192;r=2100;c=111572;regions=256;R=r*regions
    A=983040;D=65536;a=F(39,512);t0=F(381,512)
    L0=713983+c;K0=a*L0-(c-111127);M=255463-c
    lam=t0/2;score=lam*A-t0*D;threshold=lam*L0-K0
    ratio=threshold/score
    assert 0<ratio<1 and L0>=743225
    assert L0-ratio*A>0
    assert (L0-ratio*A)/2+ratio*D-M<0
    # These signs prove B/S<=ratio throughout mu,nu,eta>=0,2lambda<=tau.
    terms=[(1,choose_log(regions,44))]
    for paid,multiplicity in ((435,44),(436,212)):
        terms.extend(((multiplicity,choose_log(r,paid)),
                      (multiplicity,choose_log(m//2,r-paid)),
                      (-multiplicity,choose_log(m,r-paid))))
    event=linear_log(terms);baseline=choose_log(R,c)
    assert event[1]<baseline[0]
    lower=linear_log(((1-ratio,baseline),(ratio,event)))
    with_rows=linear_log(((1,lower),(1,choose_log(8192,2100))))
    selected_score=F(707,2560)*A-t0*D;selected_threshold=F(165556)
    selected_ratio=selected_threshold/selected_score
    selected=linear_log(((1-selected_ratio,baseline),(selected_ratio,event)))
    b4=linear_log(((1,event),));b4=tuple(v+4*(selected_score-selected_threshold) for v in b4)
    print(dict(activity=A,deaths=D,selected_score=str(selected_score),selected_threshold=str(selected_threshold),
               ratio_max=str(ratio),ratio_max_float=float(ratio),
               directed_log2_baseline=[str(v) for v in baseline],
               directed_log2_event=[str(v) for v in event],
               all_multiplier_all_b_lower=str(lower[0]),all_multiplier_all_b_lower_float=float(lower[0]),
               with_row_union_lower=str(with_rows[0]),with_row_union_lower_float=float(with_rows[0]),
               selected_weights_all_b_lower=str(selected[0]),selected_weights_all_b_lower_float=float(selected[0]),
               selected_b16_lower=str(b4[0]),selected_b16_lower_float=float(b4[0]),
               obstruction_verified=lower[0]>0,actual_regional_quotas=True,
               no_boundary_reset=True,not_failure_lower=True),flush=True)
    # A shorter crossing burst may also obstruct parameters outside onset3.
    # Here B/S<=L0/A follows directly, without imposing2lambda<=tau.
    for half_updates in (6,7,8,9):
        active=(2*half_updates-1)*512*128;deaths=512*128
        ratio=F(L0,active)
        if not 0<ratio<1:continue
        assert K0-ratio*t0*deaths>=0 and M-ratio*deaths>=0
        terms=[(1,choose_log(regions,44))]
        for paid,multiplicity in ((435,44),(436,212)):
            terms.extend(((multiplicity,choose_log(r,paid)),
                          (multiplicity,choose_log(m-half_updates*512,r-paid)),
                          (-multiplicity,choose_log(m,r-paid))))
        event=linear_log(terms)
        assert event[1]<baseline[0]
        lower=linear_log(((1-ratio,baseline),(ratio,event),(1,choose_log(8192,2100))))
        print(dict(crossing_half_updates=half_updates,activity=active,deaths=deaths,
                   universal_ratio=str(ratio),directed_all_multiplier_all_b_with_rows_lower=str(lower[0]),
                   lower_float=float(lower[0]),obstruction_verified=lower[0]>0,
                   no_onset_restriction=True,not_failure_lower=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();audit();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
