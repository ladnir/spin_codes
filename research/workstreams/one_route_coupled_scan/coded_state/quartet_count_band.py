"""FLOAT regional count-path entropy test; no finite-distance certificate.

One thread,90s watchdog,<32MiB arrays,exact regional totals, no files.
Only a bounded input-count band and payment counts435/436 are included.
"""
import itertools
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
from fractions import Fraction as F
import numpy as np
from quartet_history_screen import column_perron,coefficient_cauchy_upper,logchoose,point_mass_upper


def region_matrices(weights,initial,deviations,steps=16):
    """Rows are incoming counts,columns outgoing; extra last row starts empty."""
    n=len(deviations);lo=min(deviations);hi=max(deviations)
    state=np.zeros((hi-lo+1,n,n+1))
    for l,d in enumerate(deviations):
        state[d-lo,l,:n]=weights[:,l]
        state[d-lo,l,n]=initial[l]
    for step in range(2,steps+1):
        after=np.zeros(((hi-lo)*step+1,n,n+1))
        for l,d in enumerate(deviations):
            values=np.einsum('sai,a->si',state,weights[:,l],optimize=False)
            offset=d-lo
            after[offset:offset+state.shape[0],l,:]+=values
        state=after
    result={}
    for total in (-1,0):
        index=total-steps*lo
        result[total]=state[index].T.copy() if 0<=index<len(state) else np.zeros((n+1,n))
    return result


def tests():
    weights=np.array([[1.,2.],[3.,4.]]);initial=np.array([5.,6.]);deviations=(-1,1)
    actual=region_matrices(weights,initial,deviations,4)[0]
    direct=np.zeros((3,2))
    for start in range(3):
        for word in itertools.product(range(2),repeat=4):
            if sum(deviations[i] for i in word):continue
            value=initial[word[0]] if start==2 else weights[start,word[0]]
            for a,l in zip(word,word[1:]):value*=weights[a,l]
            direct[start,word[-1]]+=value
    assert np.array_equal(actual,direct)
    print(dict(exact_small_regional_matrix_cases=48,initial_and_terminal_indices_checked=True),flush=True)
    checked=0
    for ds in ((F(1),)*8,(F(1),F(4),F(16),F(64)),(F(1),)*3+(F(4),)*2+(F(32),)*3):
        for z in (F(1,32),F(1,4),F(1),F(4)):
            probabilities=[d*z/(1+d*z) for d in ds]
            coefficients=[F(1)]
            for p in probabilities:
                after=[F(0)]*(len(coefficients)+1)
                for k,v in enumerate(coefficients):after[k]+=v*(1-p);after[k+1]+=v*p
                coefficients=after
            variance=sum(p*(1-p) for p in probabilities)
            peak=max(coefficients)
            assert peak*peak<=F(11,28)/variance
            assert float(peak)<=point_mass_upper(float(variance),True)+1e-14
            checked+=len(coefficients)
    print(dict(exact_auxiliary_point_masses_checked=checked,
               coarse_square_bound_exact=True,refined_values_float_checked=True),flush=True)


def evaluate_band(loggamma,counts,logb,slope,tau,potentials,initial_logs=None,terminal_logs=None):
    n=len(counts);center=408;deviations=[l-center for l in counts]
    scale=float(np.max(loggamma))
    weights=np.exp2(loggamma-scale)
    if initial_logs is None:
        initial_logs=[logchoose(512,l)+l*math.log2(h[1]) for l,h in zip(counts,potentials)]
    initial=np.exp2(np.array(initial_logs)-scale)
    matrices=region_matrices(weights,initial,deviations)
    quota={-1:logchoose(6527,435)-logchoose(8192,2100),
            0:logchoose(6528,436)-logchoose(8192,2100)}
    common=max(quota.values())
    low=matrices[-1]*2**(quota[-1]-common)
    high=matrices[0]*2**(quota[0]-common)
    local_log=16*scale+common
    rewards=(slope,2*slope,3*slope,slope+tau)
    if terminal_logs is None:
        terminal_logs=[l*max(logb*r-math.log2(hi) for r,hi in zip(rewards,h[1:])) for l,h in zip(counts,potentials)]
    terminal_scale=max(terminal_logs)
    terminal=np.exp2(np.array(terminal_logs)-terminal_scale)
    # Fixed assignment,matching the preceding diagnostic's44 low regions.
    vec=low[n].copy();acc=local_log
    for region in range(1,256):
        norm=float(np.sum(vec));assert norm>0
        vec/=norm;acc+=math.log2(norm)
        vec=vec@(low[:n] if region<44 else high[:n]);acc+=local_log
    ordered=acc+math.log2(float(vec@terminal))+terminal_scale
    # Sum all choices of the44 low regions, carrying the last input count.
    with np.errstate(divide='ignore'):
        loglow=np.log2(low);loghigh=np.log2(high)
    def logmultiply(vector,matrix):
        values=vector[:,:,None]+matrix[None,:,:]
        peak=np.max(values,axis=1)
        safe=np.where(np.isfinite(peak),peak,0.)
        with np.errstate(divide='ignore'):
            return safe+np.log2(np.sum(np.exp2(values-safe[:,None,:]),axis=1))
    state=np.full((45,n),-np.inf);state[0]=loghigh[n];state[1]=loglow[n];acc=local_log
    for region in range(1,256):
        after=logmultiply(state,loghigh[:n])
        after[1:]=np.logaddexp2(after[1:],logmultiply(state[:-1],loglow[:n]))
        state=after;acc+=local_log
    final_logs=state[44]+np.array(terminal_logs)
    peak=float(np.max(final_logs));assert math.isfinite(peak)
    all_assignments=acc+peak+math.log2(float(np.sum(np.exp2(final_logs-peak))))
    shift=-logb*165556+logchoose(8192,2100)
    assert all_assignments>=ordered-1e-7
    return ordered+shift,all_assignments+shift


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();logb=4;slope=707/2560;tau=381/512;h=column_perron(logb,slope,tau)
        refined='--refined' in sys.argv
        counts=list(range(392,425)) if refined else list(range(400,417));n=len(counts)
        modes=('fixed',) if refined else ('fixed','count_adapted')
        bands=((400,416),(396,420),(392,424)) if refined else ((407,408),(407,409),(406,410),(404,412),(400,416))
        for mode in modes:
            potentials=[h for _ in counts] if mode=='fixed' else [column_perron(logb,slope,tau,a) for a in counts]
            logs=np.empty((2,n,n))
            for i,a in enumerate(counts):
                for j,l in enumerate(counts):logs[:,i,j]=coefficient_cauchy_upper(a,l,logb,slope,tau,potentials[j],True,potentials[i],point_mass='refined' if refined else 'coarse')
            print(dict(potential=mode,local_count_pairs=n*n,seconds=time.monotonic()-start),flush=True)
            for left,right in bands:
                if mode=='count_adapted' and left<404:continue  # Narrower bands already reject this candidate; avoid underflow in a pointless wider test.
                sl=slice(left-counts[0],right-counts[0]+1);band=list(range(left,right+1))
                for index,label in enumerate(('cauchy','cauchy_with_point_mass')):
                    ordered,all_assignments=evaluate_band(logs[index,sl,sl],band,logb,slope,tau,potentials[sl])
                    print(dict(input_count_band=[left,right],potential=mode,local_envelope=label,refined_point_mass=refined,
                               float_log2_ordered_regions=ordered,
                               float_log2_all_balanced_region_assignments=all_assignments,
                               regional_quotas_exact=True,total_paid=111572,
                               full_sum_outside_band_not_included=True,
                               directed_certificate=False,whole_code_distance=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
