"""FLOAT full-weight distance screen on the certified remaining count range.

Uses one global low-activity/high-turnover allowance, not one per row slice.
All numerical upper bounds still need directed verification and all-r coverage.
Serial, 120s watchdog, 128 MiB numeric-array budget, >=20 GiB free, no files.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_full_weight_kernel as full
import quartet_low_activity as low
import quartet_high_turnover as high
import quartet_quiet_regions as quiet
from quartet_regional_matrix import log_product


def midpoint_envelope(records,chunk=4):
    """Sum midpoint classes after minimizing complete-region quota bounds.

    Each record is (log half-region matrix, log quota probability). The
    same p applies to BOTH halves before the minimum. FLOAT only.
    """
    size=records[0][0].shape[0]
    result=np.full((size,size),-np.inf)
    assert size*size*8*(len(records)+2*chunk+6)<64*(1<<20)
    for lo in range(0,size,chunk):
        hi=min(size,lo+chunk);best=np.full((hi-lo,size,size),np.inf)
        for record in records:
            half,pi=record[:2]
            term=half[lo:hi,:,None]+half[None,:,:]-pi
            np.minimum(best,term,out=best)
        result[lo:hi]=np.logaddexp2.reduce(best,axis=1)
    return result


def split_quota_envelope(records,r,rho,width,degree=4096,quarters=None):
    """Keep the first-half label interval, the second's complementary range.

    Proof p may differ between halves only AFTER allocating their label
    counts. Each record includes p, and both halves use the same reference
    rho. Retains every midpoint weight. FLOAT upper envelope.
    """
    size=records[0][0].shape[0];result=np.full((size,size),-np.inf)
    assert size*size*8*(len(records)+12)<64*(1<<20)
    def envelope(lo,hi):
        best=np.full((size,size),np.inf)
        refined=[]
        for index,(half,_,p) in enumerate(records):
            slope=math.log2(rho/p)-math.log2((1-rho)/(1-p))
            endpoint=hi if slope>=0 else lo
            ratio=degree*math.log2((1-rho)/(1-p))+endpoint*slope
            np.minimum(best,half+ratio,out=best)
            if quarters is not None:refined.append((quarters[index],-ratio))
        if quarters is not None:
            best=np.minimum(best,midpoint_envelope(refined))
        return best
    bins=0
    for lo in range(max(0,r-degree),min(r,degree)+1,width):
        hi=min(min(r,degree),lo+width-1)
        left=envelope(lo,hi);right=envelope(r-hi,r-lo)
        combined,_=log_product(left,right)
        result=np.logaddexp2(result,combined);bins+=1
        if quarters is not None and size>3:
            print(dict(phase='quarter-refined quota',completed_intervals=bins,first_half_lo=lo),flush=True)
    pi=quiet.log_choose(2*degree,r)+r*math.log2(rho)+(2*degree-r)*math.log2(1-rho)
    return result-pi,bins


def bridge_weights(kernel,start,end,length):
    """Artificial positive-kernel bridge; not a law of actual code states."""
    size=len(kernel);left=[np.full(size,-np.inf)];left[0][start]=0
    for _ in range(length):left.append(np.logaddexp2.reduce(left[-1][:,None]+kernel,axis=0))
    right=[None]*(length+1);right[-1]=np.full(size,-np.inf);right[-1][end]=0
    for i in range(length-1,-1,-1):right[i]=np.logaddexp2.reduce(kernel+right[i+1][None,:],axis=1)
    normalizer=left[-1][end];weights=[]
    for i in range(length+1):
        mass=np.exp2(left[i]+right[i]-normalizer)
        assert abs(sum(mass)-1)<1e-8
        weights.append(float(mass@np.arange(size)))
    return weights


def quarter_edge(records,quarters,u,v,chunk=4):
    """One regional edge with all three interior quarter weights retained.

    One common p for all four factors, minimum before ANY interior sum.
    This is an expensive diagnostic edge, not a full operator replacement.
    """
    size=len(quarters[0]);result=-np.inf
    for lo in range(0,size,chunk):
        hi=min(size,lo+chunk);best=np.full((hi-lo,size,size),np.inf)
        for (_,pi,_),quarter in zip(records,quarters):
            left=quarter[u,:][None,:]+quarter[:,lo:hi].T
            right=quarter[lo:hi,:]+quarter[:,v][None,:]
            np.minimum(best,left[:,:,None]+right[:,None,:]-pi,out=best)
        result=np.logaddexp2(result,np.logaddexp2.reduce(best.ravel()))
    return float(result)


def tests():
    # Independent path enumeration verifies both global marker directions.
    weights=((F(1,3),F(2,7)),(F(1,5),F(3,8)))
    eta=F(3,2);zeta=F(2,3);marked=F(0);paths=[]
    for suffix in product(range(2),repeat=4):
        path=(0,)+suffix;a=sum(suffix);d=0;value=F(1)
        for u,v in zip(path,path[1:]):value*=weights[u][v];d+=int(u==1 and v==0)
        transformed=sum(path[:-1])+d+path[-1]-3
        assert transformed==a+d-3  # Free-terminal score potential telescopes.
        paths.append((a,d,value));marked+=value*eta**a*zeta**d
    vector=[F(1),F(0)]
    for _ in range(4):
        vector=[sum(vector[u]*weights[u][v]*eta**v*zeta**int(u==1 and v==0) for u in range(2)) for v in range(2)]
    assert sum(vector)==marked
    # Joint-score activity/death/coverage weights telescope only when the
    # FREE terminal contribution is also retained. Use time-varying labels
    # and signed integral coefficients for an exact independent path check.
    joint_sum=F(0);joint_good=F(0);labels=(1,0,1,0)
    for suffix in product(range(2),repeat=4):
        path=(0,)+suffix;value=F(1);c=0;d=0
        for t,(u,v) in enumerate(zip(path,path[1:])):
            value*=weights[u][v];c+=labels[t]*int(u or v);d+=int(u and not v)
        joint_score=2*sum(suffix)-d-3*c+path[-1]-2
        joint_sum+=value*F(2)**joint_score
        if joint_score>0:joint_good+=value
    vector=[F(1),F(0)]
    for label in labels:
        vector=[sum(vector[u]*weights[u][v]*F(2)**(2*v-int(u and not v)-3*label*int(u or v))
                    for u in range(2)) for v in range(2)]
    assert sum(vector[w]*F(2)**(w-2) for w in range(2))==joint_sum
    assert joint_good<=joint_sum
    for afloor in (1,2,3):
        for dcap in (0,1,2):
            exact=sum(value for a,d,value in paths if a>=afloor and d<=dcap)
            assert exact<=marked*eta**(-afloor)*zeta**(-dcap)
    # Assuming R-lambda*A >= -B, gamma^(R-lambda*A+B)>=1.
    # Integral test slope avoids floating powers in the independent check.
    cut_checks=0
    for offset in (0,1,2):
        exact=sum(value for a,d,value in paths if 2*d-a>=-offset)
        bound=sum(value*F(2)**(2*d-a+offset) for a,d,value in paths)
        assert exact<=bound;cut_checks+=1
    p=.4;theta=.175;ell=.00377;X=.23
    kernel=full.power(full.local(p,theta,ell,X,1,death_marker=.7),2)
    full.boundary_checks(kernel,p,theta,ell,X,1,2,death_marker=.7)
    joint=full.power(full.local(p,theta,ell,X,1,death_marker=.7,covered_marker=.7),2)
    full.boundary_checks(joint,p,theta,ell,X,1,2,death_marker=.7,covered_marker=.7)
    assert np.all(joint<=kernel+1e-12)
    # Every inactive position is either covered or among the free positions.
    quota_checks=0
    for size in range(1,7):
        for inactive in range(1<<size):
            for covered in range(1<<size):
                ct=(inactive&covered).bit_count()
                free=size-covered.bit_count()
                assert ct>=inactive.bit_count()-free
                quota_checks+=1
    gap=full.power(full.local(p,theta,ell,X,1,death_marker=.5,free_marker=2),2)
    full.boundary_checks(gap,p,theta,ell,X,1,2,death_marker=.5,free_marker=2)
    # Exact rational toy coefficient matrices. Condition on the full quota,
    # not separately on the two half paths.
    degree=2;size=3
    coeff=[[[F(1+i+2*j+t,17) for j in range(size)] for i in range(size)] for t in range(degree+1)]
    midpoint_checks=0;split_checks=0
    for r in (1,2,3):
        candidates=[];records=[]
        for p in (F(1,5),F(2,5),F(4,5)):
            half=[[sum(coeff[t][i][j]*p**t*(1-p)**(degree-t) for t in range(degree+1))
                   for j in range(size)] for i in range(size)]
            pi=math.comb(2*degree,r)*p**r*(1-p)**(2*degree-r)
            records.append((np.log2(np.array(half,dtype=float)),math.log2(float(pi)),float(p)))
            candidates.append([[[half[i][w]*half[w][j]/pi for j in range(size)]
                                for w in range(size)] for i in range(size)])
        got=midpoint_envelope(records)
        for i in range(size):
            for j in range(size):
                exact=sum(coeff[a][i][w]*coeff[r-a][w][j]
                          for w in range(size) for a in range(degree+1) if 0<=r-a<=degree)/math.comb(2*degree,r)
                refined=sum(min(c[i][w][j] for c in candidates) for w in range(size))
                plain=min(sum(c[i][w][j] for w in range(size)) for c in candidates)
                assert exact<=refined<=plain
                assert abs(got[i,j]-math.log2(float(refined)))<1e-12
                midpoint_checks+=1
        # Independently sum exact positive coefficients in each allocated
        # interval. The complementary interval may overcount totals, never
        # omit the exact full-region quota.
        for width in (1,2):
            split,_=split_quota_envelope(records,r,.4,width,degree)
            exact=np.zeros((size,size),dtype=object)
            for a in range(degree+1):
                if 0<=r-a<=degree:
                    exact+=np.array(coeff[a],dtype=object)@np.array(coeff[r-a],dtype=object)
            exact/=math.comb(2*degree,r)
            assert np.all(split+1e-12>=np.log2(np.array(exact,dtype=float)))
            split_checks+=size*size
    quarter_checks=0;quarter_edge_checks=0
    primitive=[np.array([[F(1+3*i+j+t,23) for j in range(size)] for i in range(size)],dtype=object) for t in (0,1)]
    exact=[np.eye(size,dtype=object)]
    for _ in range(4):
        nxt=[np.zeros((size,size),dtype=object) for _ in range(len(exact)+1)]
        for a,left in enumerate(exact):
            for b,right in enumerate(primitive):nxt[a+b]+=left@right
        exact=nxt
    for r in (1,2,3):
        records=[];quarters=[]
        for p in (F(1,5),F(2,5),F(4,5)):
            quarter=(1-p)*primitive[0]+p*primitive[1];half=quarter@quarter
            pi=math.comb(4,r)*p**r*(1-p)**(4-r)
            records.append((np.log2(np.array(half,dtype=float)),math.log2(float(pi)),float(p)))
            quarters.append(np.log2(np.array(quarter,dtype=float)))
        for width in (1,2):
            plain,_=split_quota_envelope(records,r,.4,width,2)
            refined,_=split_quota_envelope(records,r,.4,width,2,quarters)
            truth=np.log2(np.array(exact[r]/math.comb(4,r),dtype=float))
            assert np.all(refined+1e-12>=truth) and np.all(refined<=plain+1e-12)
            quarter_checks+=size*size
        midpoint=midpoint_envelope(records)
        for i in range(size):
            for j in range(size):
                edge=quarter_edge(records,quarters,i,j)
                truth=math.log2(float(exact[r][i,j]/math.comb(4,r)))
                assert truth<=edge+1e-12 and edge<=midpoint[i,j]+1e-12
                quarter_edge_checks+=1
    assert low.remaining_activity_floor(6092)==694489
    bridge=bridge_weights(np.log2(np.array(weights,dtype=float)),0,1,3)
    bridge_mass=F(0);bridge_sums=[F(0)]*4
    for middle in product((0,1),repeat=2):
        path=(0,)+middle+(1,);weight=F(1)
        for u,v in zip(path,path[1:]):weight*=weights[u][v]
        bridge_mass+=weight
        for i,w in enumerate(path):bridge_sums[i]+=w*weight
    assert all(abs(bridge[i]-float(bridge_sums[i]/bridge_mass))<1e-12 for i in range(4))
    print(dict(global_marker_path_checks=9,joint_score_terminal_paths=16,activation_reset_cells=54,
               exact_coverage_quota_checks=quota_checks,exact_midpoint_components=midpoint_checks,
               exact_split_quota_components=split_checks,exact_quarter_components=quarter_checks,
               exact_quarter_edge_components=quarter_edge_checks,
               score_telescope_paths=16,exact_bridge_marginals=4),flush=True)
    print(dict(exact_conditional_cut_checks=cut_checks,
               cut_event_requires_separate_setup_bound=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(120,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--activity',nargs='+',type=float,default=[1.,1.02,1.1,1.5])
        parser.add_argument('--death',type=float,default=1.)
        parser.add_argument('--covered',action='store_true')
        parser.add_argument('--free-gap',type=float,default=1.)
        parser.add_argument('--tau',type=float,default=0.)
        parser.add_argument('--flow',action='store_true')
        parser.add_argument('--midpoint',action='store_true')
        parser.add_argument('--split-quota',type=int,default=0)
        parser.add_argument('--quarters',action='store_true')
        parser.add_argument('--profile',action='store_true')
        parser.add_argument('--quarter-edge',action='store_true')
        parser.add_argument('--p-levels',type=int,default=9)
        parser.add_argument('--theta',type=float,default=7/40)
        parser.add_argument('--ell',type=float,default=.00377)
        parser.add_argument('--cut-marker',type=float,default=1.)
        parser.add_argument('--cut-slope',type=float,default=.24)
        parser.add_argument('--cut-offset',type=float,default=0.)
        parser.add_argument('--geometry',action='store_true')
        parser.add_argument('--geometry-marker',type=float,default=1.)
        parser.add_argument('--refine-bracket',type=int,default=-1)
        parser.add_argument('--rank-marker',type=float,default=1.)
        parser.add_argument('--markov-marker',type=float,default=1.)
        parser.add_argument('--selected-marker',type=float,default=1.)
        args=parser.parse_args();tests();low.certificate();slack=high.certificate()
        if args.tests_only:return
        n=8192-args.rows;W=943718;theta=args.theta;ell=args.ell;X=97905/432928
        assert 1<=n<8088 and 2<=args.p_levels<=9
        assert 0<theta<=1 and (1<<32)*ell>=1
        assert 0<args.death<=1 and 1<=len(args.activity)<=4 and all(x>=1 for x in args.activity)
        assert args.free_gap>=1
        assert args.cut_marker>=1 and 0<args.cut_slope<1 and args.cut_offset>=0
        assert args.geometry_marker>=1 and (args.geometry or args.geometry_marker==1)
        assert args.rank_marker>=1
        assert args.markov_marker>=1
        assert args.selected_marker>=1
        assert -1<=args.refine_bracket<8
        if args.refine_bracket>=0:
            assert args.p_levels==9 and not (args.profile or args.quarters or args.quarter_edge or args.split_quota)
        assert not (args.midpoint or args.split_quota) or len(args.activity)==1
        assert args.split_quota in (0,16,32,64,128,256,512)
        assert not args.quarters or args.split_quota>=128
        assert not args.profile or (args.midpoint and args.activity==[1.])
        assert not args.quarter_edge or (args.midpoint and args.activity==[1.])
        assert 0<=args.tau<=1
        logscore=theta*args.tau*(32+math.log2(ell))
        score=2**logscore
        afloor=low.remaining_activity_floor(n);dcap=W-128*n+slack-1
        geometry_z=0
        rank_floor=0
        markov=dict(activity=0.,deaths=0.,terminal=0.,constant=0.)
        selected=dict(activity=0.,deaths=0.,terminal=0.,constant=0.)
        if args.selected_marker>1:
            import quartet_rank_score as rank
            rank.selected_prior_certificate()
            selected={name:float(value) for name,value in rank.selected_prior_weights(n).items()}
        if args.markov_marker>1:
            import quartet_rank_score as rank
            rank.markov_certificate()
            markov={name:float(value) for name,value in rank.markov_weights(n).items()}
        if args.rank_marker>1:
            import quartet_rank_score as rank
            rank.certificate();rank_floor=rank.remaining_rank_score_floor(n)
        if args.geometry:
            from quartet_zero_neighbors import certified_row_bounds
            geometry=certified_row_bounds(args.rows)
            geometry_z=geometry['Z_upper']
            afloor=max(afloor,geometry['activity_per_active_row']*n+1)
            print(dict(certified_geometry_Z_upper=geometry_z,certified_activity_floor=afloor,
                       shared_extra_allowance='2^-60 over all unresolved row counts'),flush=True)
        free_gap_floor=256*args.rows-dcap
        regional={eta:np.full((513,513),np.inf) for eta in args.activity}
        halves=[];quarters=[];kernels=[]
        rho=args.rows/8192;odds=rho/(1-rho);next_weights=np.arange(513)[None,:]
        positions=list(range(args.p_levels))
        if args.refine_bracket>=0:positions+=[args.refine_bracket+j/9 for j in range(1,9)]
        positions.sort()
        for i in positions:
            p=1/(1+1/(odds*2**(32*theta*i/(args.p_levels-1))))
            kernel=full.power(full.local(p,theta,ell,X,score,death_marker=args.death/args.free_gap*args.cut_marker
                                        *args.markov_marker**markov['deaths']*args.selected_marker**selected['deaths'],
                                        free_marker=args.free_gap,
                                        covered_marker=(args.death if args.covered else 1)*args.cut_marker
                                        *args.geometry_marker/args.rank_marker/args.markov_marker/args.selected_marker),128)
            if args.profile:kernels.append(kernel.copy())
            pi=quiet.log_choose(8192,args.rows)+args.rows*math.log2(p)+n*math.log2(1-p)
            for eta in regional:
                direct=kernel+next_weights*(math.log2(eta)-args.cut_slope*math.log2(args.cut_marker)
                                           -math.log2(args.geometry_marker)/5+math.log2(args.rank_marker))
                direct+=next_weights*markov['activity']*math.log2(args.markov_marker)
                direct+=next_weights*selected['activity']*math.log2(args.selected_marker)
                for step in range(4):
                    if step==2 and (args.quarters or args.quarter_edge):quarters.append(direct.copy())
                    if step==3 and (args.midpoint or args.split_quota):halves.append((direct.copy(),pi,p))
                    direct,_=log_product(direct,direct)
                regional[eta]=np.minimum(regional[eta],direct-pi)
            print(dict(phase='marked regional matrices',p=p,seconds=time.monotonic()-start),flush=True)
        if args.midpoint:
            eta=args.activity[0];refined=midpoint_envelope(halves)
            assert np.max(refined-regional[eta])<1e-7
            print(dict(midpoint_improvement_min=float(np.min(regional[eta]-refined)),
                       midpoint_improvement_max=float(np.max(regional[eta]-refined)),
                       weight272_improvement=float(regional[eta][272,272]-refined[272,272]),
                       seconds=time.monotonic()-start),flush=True)
            regional[eta]=refined
            if args.quarter_edge:
                for u,v in ((0,489),(489,0)):
                    edge=quarter_edge(halves,quarters,u,v)
                    assert edge<=refined[u,v]+1e-7
                    print(dict(quarter_edge=[u,v],log2_bound=edge,
                               improvement_over_midpoint=float(refined[u,v]-edge),
                               not_a_full_operator=True,seconds=time.monotonic()-start),flush=True)
            if args.profile:
                for u,v in ((0,489),(489,0)):
                    candidates=np.array([half[u,:]+half[:,v]-pi for half,pi,p in halves])
                    scores=np.min(candidates,axis=0);w=int(np.argmax(scores));index=int(np.argmin(candidates[:,w]))
                    trace=bridge_weights(kernels[index],u,w,8)+bridge_weights(kernels[index],w,v,8)[1:]
                    print(dict(artificial_region_profile=[u,w,v],proof_p=halves[index][2],
                               midpoint_class_fraction=float(2**(scores[w]-np.logaddexp2.reduce(scores))),
                               conditional_mean_weights=[round(value,4) for value in trace],
                               actual_encoder_distribution=False),flush=True)
        if args.split_quota:
            assert (len(halves)+len(quarters)+24)*513*513*8<128*(1<<20)
            refined,bins=split_quota_envelope(halves,args.rows,rho,args.split_quota,
                                             quarters=quarters if args.quarters else None)
            eta=args.activity[0]
            print(dict(split_quota_intervals=bins,count_width=args.split_quota,
                       weight272_improvement=float(regional[eta][272,272]-refined[272,272]),
                       seconds=time.monotonic()-start),flush=True)
            regional[eta]=np.minimum(regional[eta],refined)
        base=quiet.log_choose(8192,args.rows)+theta+32*theta*(W-1048576+128*args.rows)+theta*W*math.log2(ell)
        base-=W*logscore
        for eta,matrix in regional.items():
            terminal=-np.arange(513)*(math.log2(X)+theta*math.log2(ell)-logscore)
            terminal+=np.arange(513)*markov['terminal']*math.log2(args.markov_marker)
            terminal+=np.arange(513)*selected['terminal']*math.log2(args.selected_marker)
            tail=terminal.copy()
            for _ in range(256):terminal=np.logaddexp2.reduce(matrix+terminal[None,:],axis=1)
            bound=base+terminal[0]-afloor*math.log2(eta)-dcap*math.log2(args.death)
            bound-=free_gap_floor*math.log2(args.free_gap)
            bound+=args.cut_offset*math.log2(args.cut_marker)
            bound+=geometry_z*math.log2(args.geometry_marker)/5
            bound-=rank_floor*math.log2(args.rank_marker)
            bound+=markov['constant']*math.log2(args.markov_marker)
            bound+=selected['constant']*math.log2(args.selected_marker)
            print(dict(inactive_rows=args.rows,activity_floor=afloor,remaining_score_cap=dcap,
                       eta=eta,zeta=args.death,theta=theta,ell=ell,log2_remaining_bound=float(bound),
                       zero_free_interpolation=args.tau,
                       marked_score='covered_inactive+deaths' if args.covered else 'deaths',
                       free_gap_marker=args.free_gap,free_minus_deaths_floor=free_gap_floor,
                       conditional_on_unproved_cut=args.cut_marker>1,
                       assumed_cut_slope=args.cut_slope,assumed_cut_offset=args.cut_offset,
                       certified_geometry_marker=args.geometry_marker,certified_geometry_Z=geometry_z,
                       certified_rank_marker=args.rank_marker,certified_rank_score_floor=rank_floor,
                       certified_markov_marker=args.markov_marker,
                       certified_selected_prior_marker=args.selected_marker,
                       directed=False,whole_code=False),flush=True)
            if args.flow:
                left=np.full(513,-np.inf);left[0]=0
                right=tail.copy()
                for _ in range(128):
                    left=np.logaddexp2.reduce(left[:,None]+matrix,axis=0)
                    left-=max(left)
                    before_right=right
                    right=np.logaddexp2.reduce(matrix+right[None,:],axis=1)
                    right-=max(right)
                logs=left+right;logs-=np.logaddexp2.reduce(logs)
                mass=np.exp2(logs);top=np.argsort(mass)[-8:][::-1]
                edges=left[:,None]+matrix+before_right[None,:]
                edges-=np.logaddexp2.reduce(edges.ravel());edge_mass=np.exp2(edges)
                assert np.max(abs(edge_mass.sum(axis=1)-mass))<1e-8
                dominant=np.unravel_index(np.argmax(edge_mass),edge_mass.shape)
                if args.midpoint:
                    for selected in (dominant,dominant[::-1]):
                        u,v=map(int,selected)
                        alternatives=np.array([half[u,:]+half[:,v]-pi for half,pi,p in halves])
                        scores=np.min(alternatives,axis=0);w=int(np.argmax(scores))
                        chosen=int(np.argmin(alternatives[:,w]))
                        print(dict(representative_artificial_edge=[u,v],edge_mass=float(edge_mass[u,v]),
                                   representative_midpoint=w,proof_p=halves[chosen][2],
                                   midpoint_fraction=float(2**(scores[w]-np.logaddexp2.reduce(scores))),
                                   actual_encoder_distribution=False),flush=True)
                low_indices=np.arange(128);high_indices=np.arange(385,513)
                print(dict(artificial_midpoint_flow=True,actual_encoder_distribution=False,
                           mean_boundary_weight=float(mass@np.arange(513)),
                           zero_mass=float(mass[0]),below_128=float(sum(mass[:128])),
                           above_384=float(sum(mass[385:])),
                           largest_weights=[(int(w),float(mass[w])) for w in top],
                           low_to_high=float(edge_mass[np.ix_(low_indices,high_indices)].sum()),
                           high_to_low=float(edge_mass[np.ix_(high_indices,low_indices)].sum()),
                           low_to_low=float(edge_mass[np.ix_(low_indices,low_indices)].sum()),
                           high_to_high=float(edge_mass[np.ix_(high_indices,high_indices)].sum())),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
