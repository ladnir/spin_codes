"""Fit a boundary potential to the full regional quota, not isolated updates.

FLOAT diagnostic at r2100,c111572, all counts/payments. Serial Peach,
1200s,2GiB conservative arrays,no files. No distance certificate.
"""
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
import quartet_fresh_count_screen as fresh
import quartet_fresh_regional as regional
from quartet_joint_moment_caps import geometry,edge_caps
from quartet_history_screen import logchoose


def normalize(base,scale,cells,edge_ids,labels,target_map):
    order=np.argsort(target_map);targets=target_map[order]
    starts=np.flatnonzero(np.r_[True,targets[1:]!=targets[:-1]])
    assert np.array_equal(targets[starts],np.arange(len(cells)))
    def trial(x):
        with np.errstate(divide='ignore'):logs=np.log2(base)+x*(labels-408)[None,:]
        shift=float(np.max(logs));tilted=mixed.positive_exp(logs-shift,'quota_tilt')
        matrix=np.add.reduceat(tilted[:,order],starts,axis=1)
        v=np.ones(len(cells))
        for _ in range(512):
            after=matrix@v;after/=np.max(after);after=np.maximum(after,2**-400)
            error=float(np.max(np.abs(np.log2(after/v))));v=after
            if error<1e-10:break
        rho=float(np.max((matrix@v)/v))*(1+2**-30)
        P=matrix*v[None,:]/(rho*v[:,None]);left=np.ones(len(cells))/len(cells)
        for _ in range(512):
            after=left@P;after/=sum(after);change=np.max(np.abs(after-left));left=after
            if change<1e-12:break
        edges=((left/v)@tilted)*v[target_map]/rho
        return float(edges@labels/sum(edges)),tilted,shift,v,rho,error
    lo=-2.;hi=2.
    for _ in range(14):
        x=(lo+hi)/2;result=trial(x)
        if result[0]<408:lo=x
        else:hi=x
    x=(lo+hi)/2;mean,tilted,shift,v,rho,error=trial(x)
    scale+=shift;growth=scale+math.log2(rho)
    with np.errstate(divide='ignore'):
        logs=np.log2(tilted)+np.log2(v[target_map])[None,:]-np.log2(v)[:,None]-math.log2(rho)
    P=mixed.positive_exp(logs,'quota_normalized')
    initial=np.full(len(labels),-np.inf)
    for l in np.unique(labels):
        e=edge_ids[l,l//mixed.BIN,0]
        initial[e]=logchoose(512,l)+l*math.log2(fresh.H[1])-growth+x*(l-408)+math.log2(v[target_map[e]])
    peak=max(initial);initial_v=peak+math.log2(sum(mixed.positive_exp(initial-peak,'quota_initial_sum')))
    first=mixed.positive_exp(initial-initial_v,'quota_initial')
    P=np.vstack((P,first));assert np.max(P.sum(axis=1))<=1+1e-12
    print(dict(stage='quota_normalization',marker=x,tilted_mean=mean,growth_log2=growth,
               step_vector_residual=error,step_vector_span_log2=-float(np.min(np.log2(v)))),flush=True)
    return P,x,np.log2(v),growth,initial_v


def evaluate(logrows,qs,logw,y,x):
    values=logrows+qs[:,None]*(y-x);peaks=np.max(values,axis=0)
    terms=mixed.positive_exp(values-peaks[None,:],'quota_regional_sum');sums=terms.sum(axis=0)
    output=peaks+np.log2(sums);ratios=output[:-1]-logw
    worst=int(np.argmax(ratios));rho=float(ratios[worst])
    derivative=44+256*float(qs@terms[:,worst])/float(sums[worst])
    return output,rho,derivative,worst


def fit(P,labels,target_map,cells,x,logv,growth,initial_v,iterations=8,steps=16,center=408):
    qs=np.arange(-436,1665)
    quota=logchoose(6528,436)-logchoose(8192,2100)
    quota_logs=np.array([logchoose(6528+q,436+q)-logchoose(6528,436) for q in qs])
    logw=np.zeros(len(cells));best=None;start=time.monotonic()
    terminal=mixed.terminal_logs(cells)-logv
    for iteration in range(iterations):
        # No edge enters the dummy initial cell, so its terminal value is unused.
        terminal_v=np.r_[mixed.positive_exp(logw,'quota_terminal_vector'),1.]
        rows=regional.coefficient_rows(P,labels,terminal=terminal_v,all_totals=True,steps=steps,
                                      target_map=target_map,positive_floor=2**mixed.FLOOR,center=center)
        with np.errstate(divide='ignore'):logrows=np.log2(np.array([rows[q] for q in qs]))
        logrows+=quota_logs[:,None]
        lo=x-6;hi=x
        for _ in range(8):
            if evaluate(logrows,qs,logw,lo,x)[2]<=0<=evaluate(logrows,qs,logw,hi,x)[2]:break
            width=hi-lo
            if evaluate(logrows,qs,logw,lo,x)[2]>0:lo-=width
            if evaluate(logrows,qs,logw,hi,x)[2]<0:hi+=width
        assert evaluate(logrows,qs,logw,lo,x)[2]<=0<=evaluate(logrows,qs,logw,hi,x)[2]
        for _ in range(40):
            y=(lo+hi)/2
            if evaluate(logrows,qs,logw,y,x)[2]<0:lo=y
            else:hi=y
        y=(lo+hi)/2;output,rho,_,worst=evaluate(logrows,qs,logw,y,x)
        # The initial region is evaluated directly. Only the255 later
        # regions use the maximum ratio for ordinary boundary states.
        value=initial_v+output[-1]+255*rho+4096*growth+256*quota+44*y+max(terminal-logw)
        value+=logchoose(8192,2100)-fresh.LOGB*165556
        weight=quota_logs+qs*(y-x);peak=max(weight)
        # Row sums are checked <=1+1e-12, including the rounded initial
        # row. Thirty-two floors dominate their16-step propagation.
        floor_output=mixed.FLOOR+math.log2(32)+peak+math.log2(sum(mixed.positive_exp(weight-peak,'quota_floor_sum')))
        record=dict(iteration=iteration,float_log2_upper=float(value),payment_marker=y,
                    regional_max_ratio_log2=rho,boundary_span_log2=-float(min(logw)),
                    worst_cell=cells[worst],floor_output_gap_min_bits=float(min(output)-floor_output),
                    seconds=time.monotonic()-start)
        # Lower diagnostic for this RELAXED linear transfer, not failure.
        # Remove the explicit counter floor before using coefficient rows.
        balanced={q:np.maximum(0.,rows[q]-32*2**mixed.FLOOR) for q in (-1,0)}
        if all(np.all(row>0) for row in balanced.values()):
            lower={q:np.log2(row)+quota_logs[q+436] for q,row in balanced.items()}
            lo0=float(min(lower[0][:-1]-logw));lo1=float(min(lower[-1][:-1]-logw))
            # First region has q=0. Among the255 later regions,44 have
            # q=-1 and211 have q=0. Their order is arbitrary.
            witness=initial_v+lower[0][-1]+211*lo0+44*lo1+logchoose(255,44)
            witness+=4096*growth+256*quota+44*x+min(terminal-logw)
            witness+=logchoose(8192,2100)-fresh.LOGB*165556
            record.update(float_relaxed_balanced_lower=float(witness),
                          balanced_min_ratios_log2=[lo0,lo1])
        print(record,flush=True)
        if best is None or value<best['float_log2_upper']:best=record.copy()
        after=output[:-1]-max(output[:-1]);residual=float(max(abs(after-logw)))
        if residual<1e-8:break
        logw=np.maximum((logw+after)/2,-400.);logw-=max(logw)
    counts=np.unique(labels)
    best.update(input_count_range=[int(min(counts)),int(max(counts))],all_input_counts=len(counts)==513,
                restricted_counts=None if len(counts)==513 else [int(l) for l in counts],
                all_regional_payments=True,inactive_rows=2100,covered_count=111572,
                directed_certificate=False,whole_code_distance=False,files_written=0)
    return best


def main(count_width=64,iterations=8):
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(1200,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        assert count_width in (8,16,32,64)
        intervals=[(a,min(512,a+count_width-1)) for a in range(0,513,count_width)]
        base,scale,cells,edge_ids,labels,target_map=mixed.build(intervals)
        new=edge_caps(cells,edge_ids,labels,geometry(),itertools.product((-1.,-.5,0.,.5,1.),repeat=2),progress=True)
        with np.errstate(divide='ignore'):logs=np.minimum(scale+np.log2(base),new)
        scale=float(np.max(logs));base=mixed.positive_exp(logs-scale,'quota_edge_matrix')
        P,x,logv,growth,initial_v=normalize(base,scale,cells,edge_ids,labels,target_map)
        result=fit(P,labels,target_map,cells,x,logv,growth,initial_v,iterations)
        result.update(total_seconds=time.monotonic()-start,incoming_count_width=count_width,cells=len(cells))
        print(dict(best=result),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
