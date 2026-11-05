"""FLOAT complete count-band paths with a regional, not one-step, envelope.

Peach only, one thread,90s,<768MiB working arrays,no files. No certificate.
"""
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
import numpy as np
import quartet_fresh_count_screen as fresh
from quartet_history_screen import logchoose


def build(counts):
    sizes=[len(fresh.buckets(a)) for a in counts]
    offsets=np.cumsum([0]+sizes);dimension=int(offsets[-1])
    # Conservative simultaneous dense matrices and coefficient workspaces.
    assert 6*dimension*dimension*8+4*(16*(max(counts)-min(counts))+1)*dimension*8<768*(1<<20)
    entries=[];scale=-math.inf
    for i,a in enumerate(counts):
        for j,l in enumerate(counts):
            matrix,s=fresh.transition(a,l)
            entries.append((i,j,matrix,s));scale=max(scale,s)
    matrix=np.zeros((dimension,dimension))
    for i,j,block,s in entries:
        matrix[offsets[i]:offsets[i+1],offsets[j]:offsets[j+1]]=block*2**(s-scale)
    labels=np.repeat(counts,sizes)
    return matrix,scale,labels,offsets


def tilted(base,labels,x):
    matrix=base*np.exp2(x*(labels-408))[None,:]
    v=np.ones(len(labels))
    for _ in range(512):
        after=matrix@v;after/=float(np.max(after))
        assert np.all(after>0)
        error=float(np.max(np.abs(np.log2(after/v))));v=after
        if error<1e-11:break
    # Convergence is only a proposal quality issue. The maximum row ratio
    # below remains an upper bound for every strictly positive candidate.
    right_error=error
    ratios=(matrix@v)/v;rho=float(np.max(ratios))
    probability=matrix*v[None,:]/(rho*v[:,None])
    assert np.max(np.sum(probability,axis=1))<=1+1e-12
    left=np.ones(len(labels))/len(labels)
    for _ in range(512):
        after=left@probability;after/=float(np.sum(after))
        error=float(np.max(np.abs(after-left)));left=after
        if error<1e-13:break
    return float(left@labels),v,rho,probability,right_error,error


def coefficient_rows(probability,labels,steps=16,terminal=None,all_totals=False,target_map=None,positive_floor=None,center=408):
    """Backward DP sums all count paths; row gives each incoming state.

    Each column is an edge with an exact input-count label. In the square
    case that column also identifies its destination state. target_map
    permits several exact labels to enter the same coarse destination.
    The exact total is carried before selecting a regional coefficient.
    No free terminal here.
    """
    deviations=labels-center;lo=int(min(deviations));hi=int(max(deviations))
    if target_map is None:target_map=np.arange(len(labels))
    target_map=np.asarray(target_map,dtype=int)
    dimension=probability.shape[0]
    assert len(labels)==probability.shape[1]==len(target_map)
    assert np.min(target_map)>=0 and np.max(target_map)<dimension
    state=np.ones((1,dimension)) if terminal is None else np.array(terminal)[None,:].copy()
    for step in range(1,steps+1):
        shifted=np.zeros((state.shape[0]+hi-lo,len(labels)))
        for d in np.unique(deviations):
            mask=deviations==d;offset=int(d)-lo
            shifted[offset:offset+len(state),mask]=state[:,target_map[mask]]
        state=shifted@probability.T
        if positive_floor is not None:
            # Adding coefficient mass is conservative, including at an
            # unreachable count. This optional upward floor avoids deleting
            # tiny tails in FLOAT diagnostics; it is not directed arithmetic.
            state=np.maximum(state,positive_floor)
        assert np.all(np.isfinite(state)) and np.min(state)>=0
    totals=range(steps*lo,steps*hi+1) if all_totals else (-1,0)
    return {q:state[q-steps*lo] for q in totals}


def screen(counts,fit_regions=False,all_payments=False,kernel=None):
    base,scale,labels,offsets=build(counts) if kernel is None else kernel.build(counts)
    lo=-4.;hi=4.
    for _ in range(12):
        x=(lo+hi)/2;mean,v,rho,probability,right_error,left_error=tilted(base,labels,x)
        if mean<408:lo=x
        else:hi=x
    x=(lo+hi)/2;mean,v,rho,probability,right_error,left_error=tilted(base,labels,x)
    # Add a genuine initial-state row, not a regional reset. Its potential
    # is defined by its positive tilted sum divided by the common rho.
    if kernel is None:
        initial=np.zeros(len(labels))
        for i,l in enumerate(counts):
            initial[offsets[i]+l//fresh.BIN]=2**(logchoose(512,l)+l*math.log2(fresh.H[1])-scale+x*(l-408))
    else:initial=kernel.initial_values(counts,offsets,scale,x)
    initial_v=float(initial@v)/rho
    first=initial*v/(rho*initial_v)
    augmented=np.zeros((len(labels)+1,len(labels)+1))
    augmented[:-1,:-1]=probability;augmented[-1,:-1]=first
    # There are no transitions INTO the dummy initial state, so its
    # assigned label is immaterial; use408 to keep the same counter range.
    all_labels=np.append(labels,408)
    rows=coefficient_rows(augmented,all_labels,all_totals=all_payments)
    rewards=fresh.LOGB*np.array((fresh.LAMBDA,2*fresh.LAMBDA,3*fresh.LAMBDA,fresh.LAMBDA+fresh.TAU))-np.log2(fresh.H[1:])
    terminal=[]
    if kernel is None:
        for l in counts:
            terminal.extend(max(r*rewards[0]+(l-r)*max(rewards[1:]) for r in ends) for ends in fresh.buckets(l))
    else:terminal=kernel.terminal_logs(counts)
    terminal_base=np.array(terminal)-np.log2(v)
    terminal_cost=float(np.max(terminal_base))
    growth=scale+math.log2(rho)
    if all_payments:
        # Concavity of log binomial motivates this y, but every coefficient
        # is actually evaluated; no tangent replaces a regional quota.
        logy=x-math.log2(6528/436)
        quota=logchoose(6528,436)-logchoose(8192,2100)
        qs=np.array([q for q in rows if 0<=436+q<=2100])
        with np.errstate(divide='ignore'):
            logrows=np.log2(np.array([rows[q] for q in qs]))
        logrows+=np.array([logchoose(6528+q,436+q)-logchoose(6528,436) for q in qs])[:,None]
        def payment_norm(y):
            logs=logrows+qs[:,None]*(y-x)
            peaks=np.max(logs,axis=0);weights=np.exp2(logs-peaks[None,:])
            sums=np.sum(weights,axis=0);norms=peaks+np.log2(sums)
            i=int(np.argmax(norms));meanq=float(qs@weights[:,i])/float(sums[i])
            return float(norms[i]),44+256*meanq,meanq,norms,i
        original_loggamma,_,_,_,_=payment_norm(logy)
        lo=logy-.5;hi=logy+.5
        for _ in range(48):
            mid=(lo+hi)/2;_,derivative,_,_,_=payment_norm(mid)
            if derivative<0:lo=mid
            else:hi=mid
        optimized_y=(lo+hi)/2;loggamma,_,meanq,norms,worst_index=payment_norm(optimized_y)
        payment_improvement=44*(optimized_y-logy)+256*(loggamma-original_loggamma)
        logy=optimized_y
        value=math.log2(initial_v)+4096*growth+256*quota+44*logy+256*loggamma+terminal_cost
        value+=logchoose(8192,2100)-fresh.LOGB*165556
        lower=value+256*(float(np.min(norms))-loggamma)+float(np.min(terminal_base))-terminal_cost
        balanced_lower=math.log2(initial_v)+4096*growth+44*x+float(np.min(terminal_base))
        balanced_lower+=44*(math.log2(float(np.min(rows[-1])))+logchoose(6527,435)-logchoose(8192,2100))
        balanced_lower+=212*(math.log2(float(np.min(rows[0])))+quota)
        balanced_lower+=logchoose(256,44)+logchoose(8192,2100)-fresh.LOGB*165556
        if worst_index==len(labels):worst='initial_state'
        else:
            group=int(np.searchsorted(offsets,worst_index,side='right'))-1
            if kernel is None:worst=dict(input_count=counts[group],fresh_count_bucket=fresh.buckets(counts[group])[worst_index-offsets[group]])
            else:worst=kernel.describe_label(counts[group],worst_index-offsets[group])
        return dict(input_count_band=[min(counts),max(counts)],bucket_width=fresh.BIN,
                    float_log2_all_regional_payments=value,marker_log2=x,payment_marker_log2=logy,
                    regional_weighted_log2_norm=loggamma,payment_optimization_delta=payment_improvement,
                    maximizing_row_mean_deviation=meanq,tilted_mean=mean,local_states=len(labels),
                    potential_proposal_log2_residual=right_error,
                    worst_boundary_label=worst,regional_log2_row_ratio_span=float(np.max(norms)-np.min(norms)),
                    float_log2_marked_envelope_lower=lower,float_log2_balanced_envelope_lower=balanced_lower,
                    all_count_paths_inside_band=True,all_regional_payment_allocations=True,
                    total_payments_coefficient_upper_at=111572,
                    directed_certificate=False,whole_code_distance=False)
    maxima={q:float(np.max(row)) for q,row in rows.items()}
    assert all(0<p<=1+1e-8 for p in maxima.values())
    value=math.log2(initial_v)+4096*growth+44*x+terminal_cost
    value+=44*(math.log2(maxima[-1])+logchoose(6527,435)-logchoose(8192,2100))
    value+=212*(math.log2(maxima[0])+logchoose(6528,436)-logchoose(8192,2100))
    value+=logchoose(256,44)-fresh.LOGB*165556+logchoose(8192,2100)
    if fit_regions:
        normalizers=maxima.copy();g=np.ones(len(all_labels));best=value
        constant=4096*growth+44*x
        constant+=44*(logchoose(6527,435)-logchoose(8192,2100))
        constant+=212*(logchoose(6528,436)-logchoose(8192,2100))
        constant+=logchoose(256,44)-fresh.LOGB*165556+logchoose(8192,2100)
        for iteration in range(1,13):
            after=(44/256)*rows[-1]/normalizers[-1]+(212/256)*rows[0]/normalizers[0]
            after/=float(np.max(after));assert np.all(after>0)
            residual=float(np.max(np.abs(np.log2(after/g))));g=after
            rows=coefficient_rows(augmented,all_labels,terminal=g)
            ratios={q:float(np.max(row/g)) for q,row in rows.items()}
            trial=constant+math.log2(initial_v*g[-1])+float(np.max(terminal_base-np.log2(g[:-1])))
            trial+=44*math.log2(ratios[-1])+212*math.log2(ratios[0])
            best=min(best,trial)
            print(dict(regional_fit_iteration=iteration,float_log2_family_upper=trial,
                       potential_log2_residual=residual,directed_certificate=False),flush=True)
            if residual<1e-9:break
        value=best
    return dict(input_count_band=[min(counts),max(counts)],bucket_width=fresh.BIN,
                float_log2_all_balanced_assignments=value,marker_log2=x,tilted_mean=mean,
                regional_coefficient_max=maxima,local_states=len(labels),
                potential_proposal_log2_residual=right_error,stationary_proposal_residual=left_error,
                all_count_paths_inside_band=True,regional_boundary_envelope=True,
                directed_certificate=False,whole_code_distance=False)


def tests():
    import itertools
    matrix=np.array([[.25,.5],[.75,.125]]);labels=np.array([407,409])
    checked=0
    for terminal in (np.ones(2),np.array([3.,2.])):
        actual=coefficient_rows(matrix,labels,4,terminal=terminal,all_totals=True)
        direct={q:np.zeros(2) for q in range(-4,5)}
        for initial in range(2):
            for word in itertools.product(range(2),repeat=4):
                q=sum(labels[i]-408 for i in word);previous=initial;value=1
                for j in word:value*=matrix[previous,j];previous=j
                direct[q][initial]+=value*terminal[previous];checked+=1
        assert all(np.array_equal(actual[q],direct[q]) for q in actual)
        weights={q:2.**q*math.comb(8+q,4+q) for q in actual}
        assert np.array_equal(sum(weights[q]*actual[q] for q in actual),
                              sum(weights[q]*direct[q] for q in direct))
    # A candidate need not be an eigenvector: maximum row ratio suffices.
    v=np.array([1.,3.]);rho=float(np.max((matrix@v)/v))
    normalized=matrix*v[None,:]/(rho*v[:,None])
    assert np.max(np.sum(normalized,axis=1))<=1
    assert np.all(matrix@v<=rho*v)
    print(dict(exact_backward_counter_paths=checked,weighted_regional_sums=2,
               non_eigenvector_max_ratio_checked=True),flush=True)
    matrix=np.array([[.25,.125,.5],[.125,.5,.25]])
    labels=np.array([407,408,409]);target_map=np.array([0,1,1]);checked=0
    for terminal in (np.ones(2),np.array([3.,2.])):
        actual=coefficient_rows(matrix,labels,4,terminal=terminal,all_totals=True,target_map=target_map)
        direct={q:np.zeros(2) for q in range(-4,5)}
        for initial in range(2):
            for edges in itertools.product(range(3),repeat=4):
                q=sum(labels[e]-408 for e in edges);old=initial;value=1
                for e in edges:value*=matrix[old,e];old=target_map[e]
                direct[q][initial]+=value*terminal[old];checked+=1
        assert all(np.array_equal(actual[q],direct[q]) for q in actual)
    print(dict(exact_labelled_edge_paths=checked,coarse_destinations_keep_exact_counts=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();fresh.BIN=4 if '--bin4' in sys.argv else 8
        fresh.transition.cache_clear();fresh.active_coefficients.cache_clear()
        radius=int(sys.argv[1]) if len(sys.argv)>1 else 8
        assert 1<=radius<=16
        counts=[392,407,408,424] if '--sparse-counts' in sys.argv else list(range(408-radius,409+radius))
        result=screen(counts,'--fit-regions' in sys.argv,'--all-payments' in sys.argv)
        result['allowed_input_counts']=counts
        result.update(seconds=time.monotonic()-start,files_written=0)
        print(result,flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
