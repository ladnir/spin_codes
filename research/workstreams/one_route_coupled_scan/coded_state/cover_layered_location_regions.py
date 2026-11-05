"""Integer count-region certificates for the unchanged layered construction.

No sampled point implies a region certificate. All box interiors use affine
envelopes; all active-row intervals use convexity with exact endpoint bounds.
"""
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
from math import comb
import time
import numpy as np
from layered_outer_locations import hardest,positive_rows,M,b,N
from layered_location_markers import BANK as BASE_BANK
from layered_crossover_markers import REFINED_BANK
from layered_interpolation import exponent_bounds
from layered_matching_chernoff import matching_envelope
from layered_matching_moments import factorial_envelope
from layered_two_layer_components import CONTROL_MARKERS,exact_rows
from layered_general_components import rational_operator,envelopes,mixed_dimension
from rational_log_bounds import log2_upper_scaled as log_upper_units

BANK=list(BASE_BANK)
BANK.append(('neighbor8016',24,True,F(4490228179533017,4503599627370496),[F(n,d) for n,d in (
    (1307626773173645,1152921504606846976),
    (1132855873210459,72057594037927936),
    (1960781310125813,549755813888),
    (525522860092113,137438953472),
    (1286179092988389,17179869184),
    (5712372443591307,302231454903657293676544),
    (4820638792747179,2251799813685248),
    (1524409720321031,8796093022208),
    (8812713259274457,4503599627370496),
    (702341171238089,4398046511104),
    (8175016733591901,4503599627370496),
    (5250132532720747,35184372088832),
    (7776281501023485,4503599627370496),
    (2557547986738249,17592186044416)
)]))
BANK.extend(REFINED_BANK)

# Fixed-point logarithm scale, not an exponent used for integer powering.
POWER=1<<24
SCHEDULE=((5901,6903,'s6900'),(6904,7636,'s7500'),(7637,7990,'s7936'),
          (7991,8061,'s8016'),(8062,8125,'s8064'),(8126,8192,'s8174'))


@lru_cache(maxsize=16)
def region_env(R,refined=False):
    assert refined in (False,True,'factorial')
    H,E,u=envelopes(R)
    if refined=='factorial':E=np.asarray(factorial_envelope(M,R),dtype=np.int64)
    elif refined:E=np.asarray(matching_envelope(M,R),dtype=np.int64)
    return H,E,u


@lru_cache(maxsize=16)
def local_logs(index,power=POWER):
    _,m,two,p,markers=BANK[index]
    h,rows,ends,zero=positive_rows(markers,m,p,two)
    initial=sum(log_upper_units((zero+p*h[2*c][1]*markers[5])*
                               (zero+p*h[2*c+1][1]*markers[5]),power) for c in range(2))
    return (initial,sum(log_upper_units(max(v),power) for v in rows),
            sum(log_upper_units(max(v),power) for v in ends),
            tuple(log_upper_units(1/v,power) for v in markers[:6]),
            log_upper_units(1/(1-p),power),log_upper_units((1-p)/p,power))


@lru_cache(maxsize=1)
def plain_logs():
    # All pair types are permitted; count values are irrelevant to these rows.
    profile=(192,1,1,0,0,0,0,0);m=28;markers=CONTROL_MARKERS
    h,_,ends=rational_operator(profile,markers,m,2);rows=exact_rows(profile,markers,m)
    initial=sum(log_upper_units((1+h[2*c][1]*markers[5])*(1+h[2*c+1][1]*markers[5]),POWER)
                for c in range(2))
    return initial,sum(log_upper_units(max(v),POWER) for v in rows),sum(
        log_upper_units(max(v.values()),POWER) for v in ends),tuple(log_upper_units(1/v,POWER) for v in markers[:6])


def max_linear(a,lo,hi):return sum(int(c)*(u if c>=0 else l) for c,l,u in zip(a,lo,hi))


def joint_dimension(R,lo,hi,a,rate,refined=False):
    """Upper-bound max_x a.x + rate*dbar(x), retaining shared count dependence.

    For each route-subset size choose a valid affine branch of the minimum.
    max min <= min max; this is an upper bound, not an exchange of equality.
    """
    assert rate>=0 and rate%2==0
    H,E,u=region_env(R,refined);P=(R*M-b)//2
    a=np.asarray(a,dtype=np.int64)
    D=np.array((-1,1,1,-2,0),dtype=np.int64)
    EE=np.array((-1,0,0,-1,0),dtype=np.int64)
    HH=np.array((0,1,1,-1,0),dtype=np.int64)
    base=max_linear(a,lo,hi)
    options=np.minimum.reduce((
        max_linear(a+rate*D,lo,hi)+rate*(P-128*u),
        max_linear(a+rate*EE,lo,hi)+rate*(P+E-128*u),
        max_linear(a+rate*HH,lo,hi)+rate*(H-128*u),
        base+rate*(R-128)*u))
    branch=max(base,int(options.max()))
    return min(branch,base+rate*max(0,R*M-N//2),
               max_linear(a+(rate//2)*D,lo,hi)+(rate//2)*P)


@lru_cache(maxsize=128)
def location_coeff(R,index):
    _,m,two,p,markers=BANK[index];P=(R*M-b)//2;W=400*N//1000
    ml,mu=exponent_bounds(m,POWER)
    initial,internal,terminal,mark,normal,odds=local_logs(index)
    const=POWER+mu*b+(mu-32*POWER)*P+(b//4)*(initial+terminal)
    const+=(P-b//2)//(16 if two else 4)*internal+W*mark[5]+R*M*normal
    a=tuple(v+vlog for v,vlog in zip((32*POWER-ml,-32*POWER,mu-32*POWER,32*POWER,0),mark[:5]))
    slope_s=(32*POWER-ml)*(R-128)+R*odds
    return const,a,mu,slope_s


@lru_cache(maxsize=4096)
def inverse_binomial(s):
    return log_upper_units(F(1,comb(M,s)),POWER)


def location_box(R,lo,hi,index,slo,shi,refined=False):
    assert 1<=slo<=shi<=M
    const,a,rate,slope_s=location_coeff(R,index)
    fixed=const+joint_dimension(R,lo,hi,a,rate,refined)
    endpoints=[fixed+s*slope_s+(R-1)*inverse_binomial(s) for s in {slo,shi}]
    return -((-max(endpoints))//POWER)


def plain_box(R,lo,hi,cap=5900,refined=False):
    # Union over ALL active sets of cardinality <=cap costs at most 2^M.
    m=28;P=(R*M-b)//2;initial,internal,terminal,mark=plain_logs()
    const=POWER*(1+m*b+(m-32)*P+(32-m)*(R-128)*cap+M)
    const+=(b//4)*(initial+terminal)+(P-b//2)//16*internal+(400*N//1000)*mark[5]
    a=tuple(POWER*v+vlog for v,vlog in zip((32-m,-32,m-32,32,0),mark[:5]))
    return -((-(const+joint_dimension(R,lo,hi,a,m*POWER,refined)))//POWER)


def count_bits(lo,hi):
    n=1
    for a,c in zip(lo,hi):n*=c-a+1
    return (n-1).bit_length()


def cover_box(R,lo,hi):
    assert 129<=R<=256 and all(0<=l<=h<=N for l,h in zip(lo,hi)) and len(lo)==len(hi)==5
    pieces=[(1,5900,'plain',plain_box(R,lo,hi))];cursor=5901
    for left,right,name in SCHEDULE:
        assert left==cursor
        index=next(i for i,v in enumerate(BANK) if v[0]==name)
        pieces.append((left,right,name,location_box(R,lo,hi,index,left,right)));cursor=right+1
    assert cursor==M+1
    per_class=max(v[3] for v in pieces)+13
    return dict(per_class=per_class,box_union=per_class+count_bits(lo,hi),pieces=pieces)


def adaptive_box(R,lo,hi,step=32,refined=False):
    """Pick proof weights per active-row interval, then replay exactly."""
    assert 129<=R<=256 and all(0<=l<=h<=N for l,h in zip(lo,hi)) and len(lo)==len(hi)==5
    assert isinstance(step,int) and step>0
    pieces=[(1,5900,'plain',plain_box(R,lo,hi,refined=refined))]
    constants=[];slopes=[]
    for i in range(len(BANK)):
        c,a,rate,slope=location_coeff(R,i)
        constants.append(c+joint_dimension(R,lo,hi,a,rate,refined));slopes.append(slope)
    constants=np.asarray(constants,dtype=np.int64);slopes=np.asarray(slopes,dtype=np.int64)
    for left in range(5901,M+1,step):
        right=min(M,left+step-1)
        # Integer endpoint majorants. Choosing the best bank is exact too.
        at_left=constants+slopes*left+(R-1)*inverse_binomial(left)
        at_right=constants+slopes*right+(R-1)*inverse_binomial(right)
        scores=np.maximum(at_left,at_right);i=int(scores.argmin())
        value=-((-int(scores[i]))//POWER)
        if len(pieces)>1 and pieces[-1][2]==BANK[i][0]:
            a,c,name,old=pieces[-1];pieces[-1]=(a,right,name,max(old,value))
        else:pieces.append((left,right,BANK[i][0],value))
    per_class=max(v[3] for v in pieces)+13
    return dict(per_class=per_class,box_union=per_class+count_bits(lo,hi),pieces=pieces)


def tube(R,zmax,radius,budget=256):
    """Certify a full integer tube by subdivision; never discard failed cells."""
    assert isinstance(zmax,int) and zmax>=0 and isinstance(radius,int) and radius>=0
    assert isinstance(budget,int) and budget>0
    center=hardest()[1:6];lo_tail=tuple(v-radius for v in center[1:])
    hi_tail=tuple(v+radius for v in center[1:]);pending=[(0,zmax)];leaves=[];failed=[];tested=0
    while pending:
        left,right=pending.pop();lo=(left,)+lo_tail;hi=(right,)+hi_tail
        result=cover_box(R,lo,hi);tested+=1
        if result['box_union']<-100:leaves.append((left,right,result['per_class'],result['box_union']))
        elif left<right and tested+len(pending)+2<=budget:
            mid=(left+right)//2;pending.extend(((mid+1,right),(left,mid)))
        else:failed.append((left,right,result['box_union']))
    leaves.sort();failed.sort();all_intervals=sorted([(a,c) for a,c,_,_ in leaves]+[(a,c) for a,c,_ in failed])
    cursor=0
    for a,c in all_intervals:assert a==cursor;cursor=c+1
    assert cursor==zmax+1
    union=max((v[3] for v in leaves),default=0)+(max(1,len(leaves))-1).bit_length()
    failures=[]
    for left,right,value in failed:
        if failures and failures[-1][1]+1==left:
            a,c,v=failures[-1];failures[-1]=(a,right,max(v,value))
        else:failures.append((left,right,value))
    return dict(R=R,z=(0,zmax),side_radius=radius,certified_boxes=len(leaves),tested=tested,
                failed=failures,failed_cells=len(failed),complete=not failed,log2_union_upper=union if not failed else None,
                leaves=leaves)


def subdivide_box(R,lo,hi,budget=256,step=8,refined=False):
    """Gapless multidimensional cover; keep every failed leaf explicitly."""
    assert len(lo)==len(hi)==5 and all(0<=a<=c<=N for a,c in zip(lo,hi))
    assert isinstance(budget,int) and budget>0
    pending=[(tuple(lo),tuple(hi))];accepted=[];failed=[];tested=0
    def volume(a,c):
        v=1
        for l,h in zip(a,c):v*=h-l+1
        return v
    while pending:
        a,c=pending.pop();result=adaptive_box(R,a,c,step,refined);tested+=1
        if result['box_union']<-100:
            accepted.append((a,c,result['box_union']))
        elif a!=c and tested+len(pending)+2<=budget:
            axis=max(range(5),key=lambda j:c[j]-a[j]);mid=(a[axis]+c[axis])//2
            left=list(c);left[axis]=mid;right=list(a);right[axis]=mid+1
            pending.extend(((tuple(right),c),(a,tuple(left))))
        else:failed.append((a,c,result['box_union']))
    assert tested<=budget
    assert sum(volume(a,c) for a,c,_ in accepted+failed)==volume(lo,hi)
    # The binary split itself proves disjoint coverage. Also check the final
    # leaves independently, so an accidental overlap cannot hide a hole.
    leaves=accepted+failed
    for i,(a,c,_) in enumerate(leaves):
        for aa,cc,_ in leaves[:i]:
            assert any(h<l or hh<ll for ll,h,l,hh in zip(a,c,aa,cc))
    bound=max((v for _,_,v in accepted),default=0)+(max(1,len(accepted))-1).bit_length()
    return dict(R=R,lo=lo,hi=hi,refined_route=refined,tested=tested,accepted_boxes=len(accepted),
        failed_boxes=len(failed),complete=not failed,
        log2_union_upper=bound if not failed else None,accepted=accepted,failed=failed)


def verify_crossover():
    """Replay selected neighbors AND full regions; report them separately."""
    began=time.monotonic();c=hardest()[1:6];neighbors=[]
    for ds in product((-128,0,128),repeat=5):
        x=tuple(a+d for a,d in zip(c,ds))
        neighbors.append((adaptive_box(192,x,x,step=1)['per_class'],ds))
    assert len(neighbors)==243 and max(v for v,_ in neighbors)<-100
    lo=(0,)+tuple(v-128 for v in c[1:]);hi=tuple(v+128 for v in c)
    regions=[subdivide_box(192,lo,hi),subdivide_box(193,lo,hi,refined=True),
        subdivide_box(192,(0,)+tuple(v-256 for v in c[1:]),
                      (50000,)+tuple(v+256 for v in c[1:]),refined=True)]
    assert all(v['complete'] and v['log2_union_upper']<-100 for v in regions)
    dims=[mixed_dimension(193,(193*M-b)//2,*c[:4],region_env(193,v)) for v in (False,True)]
    frontier=[(R,adaptive_box(R,c,c,step=1,refined=True)['per_class']) for R in (193,194,195)]
    print(dict(selected_neighbors=len(neighbors),worst_neighbor=max(neighbors),
        neighbor_union_upper=max(v for v,_ in neighbors)+8,
        regions=[{k:x for k,x in v.items() if k not in ('accepted','failed')} for v in regions],
        R193_dimensions=dims,refined_route_frontier=frontier,
        elapsed_seconds=time.monotonic()-began,whole_code_certificate=False))


def verify_suffix_progress():
    """Full count-box replays under the actual matching-moment envelope."""
    began=time.monotonic();c=hardest()[1:6];regions=[]
    for R,z,radius in ((193,51000,512),(194,50000,256),(195,49264,128)):
        lo=(0,)+tuple(v-radius for v in c[1:]);hi=(z,)+tuple(v+radius for v in c[1:])
        result=subdivide_box(R,lo,hi,refined='factorial')
        assert result['complete'] and result['log2_union_upper']<-100
        regions.append({k:v for k,v in result.items() if k not in ('accepted','failed')})
    dims=[(R,[mixed_dimension(R,(R*M-b)//2,*c[:4],region_env(R,v))
              for v in (True,'factorial')]) for R in (193,194,195)]
    frontier=[(R,adaptive_box(R,c,c,step=1,refined='factorial')['per_class']) for R in range(193,199)]
    print(dict(regions=regions,dimensions=dims,central_profiles=frontier,
        elapsed_seconds=time.monotonic()-began,whole_code_certificate=False))


def run(neighbors=False,adaptive=False,tubes=False,active_step=32):
    began=time.monotonic();p=hardest();center=p[1:6];values=[]
    if tubes:
        results=[tube(192,49200,64),tube(192,48000,128),tube(190,49168,128),tube(191,49168,128)]
        assert all(v['complete'] and v['log2_union_upper']<-100 for v in results)
        adjacent=[]
        for R in (190,191,193,194):
            adjacent.append((R,cover_box(R,center,center)['per_class']))
        print(dict(tubes=results,adjacent_suffix_points=adjacent,elapsed_seconds=time.monotonic()-began,
                   whole_code_certificate=False));return
    if adaptive:
        tested=[]
        for delta in ((0,0,0,0,0),(128,-128,-128,128,-128)):
            x=tuple(v+d for v,d in zip(center,delta));tested.append((delta,adaptive_box(192,x,x,active_step)))
        sizes=[]
        for radius in (32,64,128,256):
            lo=tuple(v-radius for v in center);hi=tuple(v+radius for v in center)
            result=adaptive_box(192,lo,hi,active_step);sizes.append((radius,result['box_union']))
        print(dict(points=tested,boxes=sizes,elapsed_seconds=time.monotonic()-began,whole_code_certificate=False))
        return
    if neighbors:
        tested=[]
        for delta in product((-128,0,128),repeat=5):
            x=tuple(v+d for v,d in zip(center,delta));value=cover_box(192,x,x)['per_class']
            tested.append((value,delta))
        zslice=[]
        for z in range(0,98305,4096):
            x=(z,)+center[1:];zslice.append((z,cover_box(192,x,x)['per_class']))
        print(dict(selected_neighbors=len(tested),passed=sum(v<-100 for v,_ in tested),
                   worst=max(tested),best=min(tested),z_slice=zslice,
                   elapsed_seconds=time.monotonic()-began,whole_code_certificate=False))
        return
    for radius in (0,16,32,64,128,256,512,1024):
        lo=tuple(v-radius for v in center);hi=tuple(v+radius for v in center)
        result=cover_box(192,lo,hi)
        values.append((radius,result['per_class'],result['box_union']))
    print(dict(center=center,box_radius_screen=values,elapsed_seconds=time.monotonic()-began,
               whole_code_certificate=False))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--neighbors',action='store_true');p.add_argument('--adaptive',action='store_true');p.add_argument('--tubes',action='store_true');p.add_argument('--crossover',action='store_true');p.add_argument('--suffix-progress',action='store_true');p.add_argument('--active-step',type=int,default=32);args=p.parse_args()
    if args.suffix_progress:verify_suffix_progress()
    elif args.crossover:verify_crossover()
    else:run(args.neighbors,args.adaptive,args.tubes,args.active_step)
