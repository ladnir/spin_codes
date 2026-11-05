"""Demand-dependent allocation of target-gate capacity to two timelines.

Construction unchanged. Exact integer timeline cuts and monotone local
allocation; selected-window diagnostic, not a setup or distance proof.
Standalone standard library, serial <=30s, no data files.
"""
import random
import shutil
import tempfile
import time
from itertools import product
from fractions import Fraction as F


def line_cut(service,supply,edge):
    f0=f1=0
    for capacity in service:
        f0,f1=min(f0,f1+edge)+supply,min(f1,f0+edge)+capacity
    return min(f0,f1)


def threshold(service,index,supply,edge):
    # Minimum costs outside this node, conditional on its two cut sides.
    l0=l1=0
    for j in range(index):
        l0,l1=min(l0,l1+edge)+supply,min(l1,l0+edge)+service[j]
    if index:l0,l1=min(l0,l1+edge),min(l1,l0+edge)
    r0=r1=0
    for j in range(len(service)-1,index,-1):
        r0,r1=min(supply+r0,service[j]+edge+r1),min(supply+edge+r0,service[j]+r1)
    return l0+supply+r0-l1-r1


def allocate(old,total,a,b):
    a=max(0,a);b=max(0,b)
    if a+b<=total:low,high=a,total-b
    else:low,high=max(0,total-b),min(total,a)
    new=max(low,min(high,old))
    gain=min(new,a)+min(total-new,b)-min(old,a)-min(total-old,b)
    assert gain>=0 and 0<=new<=total
    return new,gain


def causal_pair(a,b,targets):
    """Scale4, unit supply, temporal capacity2; five-state causal rule."""
    total=4*targets
    left=max(0,min(total,(total+b-a)//2));right=total-left
    gain=max(-a,0)+min(1,a+left)+max(-b,0)+min(1,b+right)
    na=max(-2,min(2,a+left-1));nb=max(-2,min(2,b+right-1))
    loss=(int(a==-2)+int(b==-2)) if not targets else 0
    assert gain==2+max(-a,0)+max(-b,0)-max(-na,0)-max(-nb,0)-loss
    return na,nb,gain,loss


def mean_field_causal(p,steps=4096):
    # Independent partner-state approximation; not the finite ensemble law.
    law=[0.,0.,1.,0.,0.];flow=0.;losses=0.
    target_law=((0,(1-p)**2),(1,2*p*(1-p)),(2,p*p))
    for step in range(steps):
        new=[0.]*5
        choices=target_law if step<steps-1 else ((0,1.),)
        for a in range(-2,3):
            for b in range(-2,3):
                for targets,chance in choices:
                    weight=law[a+2]*law[b+2]*chance
                    na,nb,gain,loss=causal_pair(a,b,targets)
                    new[na+2]+=weight/2;new[nb+2]+=weight/2
                    flow+=weight*gain/2;losses+=weight*loss/2
        # Avoid amplifying binary floating point mass drift under the quadratic map.
        mass=sum(new);law=[x/mass for x in new]
    return dict(mean_field_only=True,steps=steps,state_law=law,
                mean_flow_per_coordinate=flow/4,mean_loss_per_step=losses/steps,
                extrapolated_full_dual=(flow/4-steps/8)*512,
                finite_population_or_setup_bound=False)


def fixed_pair_stationary():
    p=F(1974,8192);target_law=((0,(1-p)**2),(1,2*p*(1-p)),(2,p*p))
    matrix=[[F(0) for _ in range(5)] for _ in range(5)]
    for a in range(-2,3):
        for targets,chance in target_law:
            na,nb,_,_=causal_pair(a,a,targets)
            assert na==nb
            matrix[a+2][na+2]+=chance
    equations=[[matrix[j][i]-int(i==j) for j in range(5)]+[F(0)] for i in range(4)]
    equations.append([F(1)]*5+[F(1)])
    for j in range(5):
        pivot=next(i for i in range(j,5) if equations[i][j])
        equations[j],equations[pivot]=equations[pivot],equations[j]
        factor=equations[j][j];equations[j]=[x/factor for x in equations[j]]
        for i in range(5):
            if i!=j:
                factor=equations[i][j]
                equations[i]=[x-factor*y for x,y in zip(equations[i],equations[j])]
    law=[equations[i][-1] for i in range(5)]
    assert sum(law)==1 and all(x>=0 for x in law)
    assert all(sum(law[j]*matrix[j][i] for j in range(5))==law[i] for i in range(5))
    loss=law[0]*(1-p)**2
    return dict(fixed_pair_iid_control=True,exact_rational_stationarity_checked=True,
                stationary_law=[float(x) for x in law],loss_per_coordinate=float(loss),
                critical_loss_fraction=371901/(4095*512),
                construction_counterexample=False)


def improve(service,pairs,gate_capacity,supply,edge,rounds=4):
    n=len(service);L=len(service[0]);flow=sum(line_cut(s,supply,edge) for s in service)
    records=[]
    for iteration in range(rounds):
        gain=0;changes=0
        order=range(L-1) if iteration%2==0 else range(L-2,-1,-1)
        for t in order:
            demands=[threshold(s,t,supply,edge) for s in service]
            for (a,b),total in zip(pairs[t],gate_capacity[t]):
                assert service[a][t]+service[b][t]==total
                new,added=allocate(service[a][t],total,demands[a],demands[b])
                changes+=new!=service[a][t];gain+=added
                service[a][t]=new;service[b][t]=total-new
        checked=sum(line_cut(s,supply,edge) for s in service)
        assert checked==flow+gain,'local gain must equal the exact timeline-cut gain'
        flow=checked;records.append(dict(round=iteration,gain=gain,changes=changes,flow=flow))
        if not changes:break
    return flow,records


def self_test():
    thresholds=allocations=0
    for capacities in product((0,2,4),repeat=4):
        for index in range(4):
            demand=threshold(capacities,index,1,2)
            baseline=None
            for value in range(9):
                trial=list(capacities);trial[index]=value
                cut=line_cut(trial,1,2)
                offset=cut-min(value,demand)
                if baseline is None:baseline=offset
                assert offset==baseline
                thresholds+=1
    for total in (0,4,8):
        for a,b in product(range(-3,6),repeat=2):
            optimum=max(min(j,max(a,0))+min(total-j,max(b,0)) for j in range(total+1))
            for old in range(total+1):
                new,gain=allocate(old,total,a,b)
                assert min(new,max(a,0))+min(total-new,max(b,0))==optimum
                allocations+=1
    # Exhaustive cuts independently check the DP itself.
    cuts=0
    for capacities in product((0,2,4),repeat=4):
        brute=min(sum((1-((mask>>j)&1))+capacities[j]*((mask>>j)&1) for j in range(4))
                  +2*sum(((mask>>j)^(mask>>(j-1)))&1 for j in range(1,4))
                  for mask in range(16))
        assert line_cut(capacities,1,2)==brute;cuts+=1
    stability=0
    for a,b,targets in product(range(-2,3),range(-2,3),range(3)):
        low=causal_pair(a,b,targets)
        for which in (0,1):
            values=[a,b]
            if values[which]==2:continue
            values[which]+=1;high=causal_pair(*values,targets)
            assert high[0]>=low[0] and high[1]>=low[1]
            assert 0<=low[3]-high[3]
            assert high[0]-low[0]+high[1]-low[1]+low[3]-high[3]<=1
            stability+=1
        if targets<2:
            high=causal_pair(a,b,targets+1)
            assert high[0]>=low[0] and high[1]>=low[1]
            assert high[0]-low[0]+high[1]-low[1]+low[3]-high[3]<=4
            stability+=1
    return dict(threshold_checks=thresholds,allocation_checks=allocations,exact_line_cuts=cuts,
                causal_stability_checks=stability)


def variable_weight_flow_tests():
    """All support histories for small fixed geometries, including births/resets.

    The cut bound depends on aggregate activity and endpoint weights, not
    on constant weight at each update. This is a deterministic geometric
    interface; no setup-probability or field-distance claim is tested here.
    """
    n=4;length=3;checks=0;setups=0;nonconstant=0
    matchings=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
    for first,second in product(range(3),repeat=2):
        layers=(matchings[first],matchings[second])
        # Include empty/full tags and asymmetric labels. Pair capacities
        # and per-coordinate coverage are computed independently below.
        for targets in ((0,0,0),(15,15,15),(5,3,10),(10,12,5)):
            aux=[0]*n;service=[[0]*length for _ in range(n)]
            random_loss=0;gain_total=0
            for t,pairs in enumerate(layers):
                for a,b in pairs:
                    j=((targets[t+1]>>a)&1)+((targets[t+1]>>b)&1)
                    left=max(0,min(4*j,(4*j+aux[b]-aux[a])//2))
                    service[a][t]=left;service[b][t]=4*j-left
                    na,nb,gain,loss=causal_pair(aux[a],aux[b],j)
                    aux[a],aux[b]=na,nb;gain_total+=gain;random_loss+=loss
            terminal_loss=0
            for a,b in matchings[0]:
                na,nb,gain,loss=causal_pair(aux[a],aux[b],0)
                aux[a],aux[b]=na,nb;gain_total+=gain;terminal_loss+=loss
            negative=sum(max(-x,0) for x in aux)
            flow=sum(line_cut(row,1,2) for row in service)
            assert flow==gain_total==n*length-random_loss-terminal_loss-negative
            assert terminal_loss+negative<=3*n;checks+=2;setups+=1
            for masks in product(range(1<<n),repeat=length):
                weights=[mask.bit_count() for mask in masks];activity=sum(weights)
                deaths=sum((a&~b).bit_count() for a,b in zip(masks,masks[1:]))
                variation=sum((a^b).bit_count() for a,b in zip(masks,masks[1:]))
                assert 2*deaths==variation+weights[0]-weights[-1]
                coverage=(targets[0]&masks[0]).bit_count()
                for t,pairs in enumerate(layers):
                    live=0
                    for a,b in pairs:
                        pair=(1<<a)|(1<<b)
                        if masks[t]&pair:live|=pair
                    coverage+=(targets[t+1]&(live|masks[t+1])).bit_count()
                service_cost=sum(service[i][t]*((masks[t]>>i)&1) for i in range(n) for t in range(length))
                assert service_cost<=4*coverage
                assert flow<=n*length-activity+service_cost+2*variation
                adjusted_activity=activity+2*(weights[0]-weights[-1])
                assert 4*(coverage+deaths)>=flow-n*length+adjusted_activity
                assert 4*(coverage+deaths)>=adjusted_activity-random_loss-3*n
                nonconstant+=int(len(set(weights))>1);checks+=5
    assert nonconstant>0
    # Same threshold as the restricted half-full calculation when the
    # adjusted activity is at least half of all coordinate-updates.
    activity=512*4096//2;required=168785
    loss_limit=activity-3*512-4*required
    assert loss_limit==371900
    assert loss_limit-2*512==370876;checks+=2
    return dict(variable_weight_checks=checks,setups=setups,nonconstant_histories=nonconstant,
                same_adjusted_activity_loss_limit=loss_limit,
                aggregate_half_activity_worst_endpoint_loss_limit=loss_limit-2*512,
                whole_code_certificate=False)


def screen(regional,seed=20261028):
    rng=random.Random(seed);n=512;L=32;scale=4;supply=1;edge=2;p=1974/8192
    if regional:
        flat=[]
        for _ in range(2):
            region=[True]*1974+[False]*(8192-1974);rng.shuffle(region);flat.extend(region)
        target=[flat[t*n:(t+1)*n] for t in range(L)]
    else:target=[[rng.random()<p for _ in range(n)] for _ in range(L)]
    pairs=[];capacities=[];service=[[0]*L for _ in range(n)]
    for t in range(L-1):
        order=list(range(n));rng.shuffle(order)
        layer=list(zip(order[::2],order[1::2]));pairs.append(layer);caps=[]
        for a,b in layer:
            total=scale*(int(target[t+1][a])+int(target[t+1][b]));caps.append(total)
            service[a][t]=service[b][t]=total//2
        capacities.append(caps)
    equal=sum(line_cut(s,supply,edge) for s in service)
    final,records=improve(service,pairs,capacities,supply,edge)
    def lower(value):return (value-supply*L*n/2)/scale
    return dict(regional_fixed_count=regional,seed=seed,window=L,
                equal_cut_lower=lower(equal),demand_cut_lower=lower(final),rounds=records,
                extrapolated_equal_cost=lower(equal)*4096/L,
                extrapolated_demand_cost=lower(final)*4096/L,
                actual_full_geometry_covered=False,whole_code_certificate=False)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--causal',action='store_true')
    parser.add_argument('--tests-only',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();print(self_test(),flush=True)
    if args.tests_only:
        print(variable_weight_flow_tests(),flush=True)
        assert time.monotonic()-began<25
        print(dict(seconds=time.monotonic()-began));return
    if args.causal:
        print(mean_field_causal(1974/8192),flush=True)
        print(fixed_pair_stationary(),flush=True)
        print(dict(seconds=time.monotonic()-began,no_setup_tail_bound=True));return
    for regional in (False,True):
        print(screen(regional),flush=True)
        assert time.monotonic()-began<25
    print(dict(seconds=time.monotonic()-began,no_setup_tail_bound=True))


if __name__=='__main__':main()
