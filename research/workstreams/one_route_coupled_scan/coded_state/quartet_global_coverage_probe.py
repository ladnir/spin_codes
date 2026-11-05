"""Global rank diagnostics; floating evaluations are NOT certificates.

No encoder changes. Exact checks use small widths and rational arithmetic.
The numerical route census retains every regional coverage quota. Run on
Peach, file-free, with one thread, <=180s and <=2 GiB address space.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, lgamma, log, log2, prod
import os
import resource
import shutil
import sys
import threading
import time

N=1<<21
M=8192
K=1<<20
W=943718


def poly_mul(left,right):
    out={}
    for (u,v),a in left.items():
        for (uu,vv),b in right.items():
            key=(u+uu,v+vv)
            out[key]=out.get(key,F(0))+a*b
    return out


def transfer_row(width,u,field_marker,output_marker,coverage_marker,survivor_marker=None):
    """Exact marked all-destination row, averaged over uniform quartets.

    Weight is a^(j-deaths) z^j x^|V|, not the whole routing bound.
    Here a is a positive field moment, 0<z<=1, x>0.
    """
    assert width%4==0
    a=F(field_marker); z=F(output_marker); x=F(coverage_marker); m=a*z
    survivor=1+m if survivor_marker is None else F(survivor_marker)
    local={}
    for incoming in range(5):
        for deaths in range(incoming+1):
            survivors=incoming-deaths
            for births in range(5-incoming):
                outgoing=survivors+births
                covered=4 if incoming else births
                weight=(comb(4,incoming)*comb(incoming,deaths)
                        *comb(4-incoming,births)*z**deaths*m**births
                        *survivor**survivors*x**covered)
                key=(incoming,outgoing)
                local[key]=local.get(key,F(0))+weight
    total={(0,0):F(1)}
    for _ in range(width//4): total=poly_mul(total,local)
    return {v:total.get((u,v),F(0))/comb(width,u) for v in range(width+1)}


def partitions(items):
    if not items:
        yield ()
        return
    first=items[0]
    for tail in combinations(items[1:],3):
        group=(first,)+tail
        rest=tuple(i for i in items if i not in group)
        for other in partitions(rest): yield (group,)+other


def exact_transfer_checks(group_outputs=False):
    checks=0
    if group_outputs:
        # q=16, theta=1/2: a=4,z=3/16,m=3/4,c=sqrt(1+m^2)=5/4.
        a=F(4);z=F(3,16);x=F(5,6);survivor=F(5,4)
        assert survivor**2==1+(a*z)**2
    else:
        a=F(3);z=F(2,7);x=F(5,6);survivor=1+a*z
    for width in (4,8):
        parts=list(partitions(tuple(range(width))))
        for u in range(width+1):
            incoming=(1<<u)-1
            exact={v:F(0) for v in range(width+1)}
            for part in parts:
                live=0
                for group in part:
                    bits=sum(1<<i for i in group)
                    if incoming&bits: live|=bits
                for outgoing in range(1<<width):
                    births=(outgoing&~incoming).bit_count()
                    deaths=(incoming&~outgoing).bit_count()
                    survivors=(incoming&outgoing).bit_count()
                    covered=(live|outgoing).bit_count()
                    # Sum all optional nonzero output positions on survivors.
                    weight=(a*z)**births*z**deaths*survivor**survivors*x**covered
                    exact[outgoing.bit_count()]+=weight/len(parts)
            predicted=transfer_row(width,u,a,z,x,survivor)
            for v in range(width+1):
                assert exact[v]==predicted[v]
                checks+=1
    print(dict(exact_marked_all_destination_checks=checks,outputs_grouped_before_moment=group_outputs),flush=True)


def rational_coverage_certificate():
    # c=2^(1/4)<6/5. Potential (1,1/4), eigenvalue upper5/4.
    assert F(6,5)**4>2
    assert F(1,256)+F(6,5)*F(1,4)<F(5,4)*F(1,4)
    assert F(1)+F(1,4)==F(5,4)
    # rho=16/7,p=7/8,t=256 balances the two positive MGF terms.
    p=F(7,8);rho=F(16,7);t=256
    assert 1-p+p*rho==F(17,8)
    assert t**128*(1-p+p*rho/t)**256==F(17,8)**256
    # No logarithmic rounding is used in the reported certificate.
    assert M+1<1<<14
    assert 85**128*7**112<1<<1135
    exponent=-7+8*(W-K)+M+1024+14*256+47*N//128
    assert exponent==-56023
    print(dict(rational_high_coverage_threshold=7*N//8,
               all_history_failure_strictly_below_2_to=exponent,
               full_minimum_distance_certificate=False),flush=True)
    # A stronger rational specialization, without a numerical optimizer.
    c=F(609,512);h=F(107,512);g=F(619,512)
    assert c**4>2
    assert 1+h==g and F(1,256)+c*h<g*h
    assert 1/h<8  # endpoint factor <=2^(3*512)
    p=F(53,64);rho=F(176,53)
    assert 1-p+p*rho==F(187,64)
    assert 256**128*(1-p+p*rho/256)**256==F(187,64)**256
    assert 115753**64*53**53<(1<<985)*176**53
    exponent=-7+8*(W-K)+M+1536+14*256+25*N//64
    assert exponent==-6359
    print(dict(rational_improved_coverage_threshold=53*N//64,
               all_history_failure_strictly_below_2_to=exponent,
               full_minimum_distance_certificate=False),flush=True)


def exact_fringe_checks():
    # Four regions, three rows, cap two. Pin all births in region zero.
    # Conditional remaining coverage in that region is uniform as well.
    rows=3;regions=4;cap=2;a=F(3);checks=0
    for coverage in product(range(rows+1),repeat=regions):
        for births in range(1,coverage[0]+1):
            required=set(range(births))
            first=[required|set(tail) for tail in combinations(range(births,rows),coverage[0]-births)]
            choices=[first]+[[set(x) for x in combinations(range(rows),v)] for v in coverage[1:]]
            numerator=F(0);denominator=1
            for sets in choices:denominator*=len(sets)
            for selected in product(*choices):
                degrees=[sum(r in s for s in selected) for r in range(rows)]
                if all(degrees[r]>=cap+1 for r in required):
                    numerator+=a**(-sum(min(cap,v) for v in degrees))
            exact=numerator/denominator
            for size in (1,2):
                for tail in combinations(range(1,regions),size-1):
                    core={0,*tail};needed=cap+1-size
                    probabilities=[F(coverage[r],rows) for r in range(regions) if r not in core]
                    elementary=sum((prod(terms,start=F(1))
                                    for terms in combinations(probabilities,needed)),F(0))
                    upper=a**(-sum(coverage[r] for r in core))*min(F(1),elementary**births)
                    assert exact<=upper
                    checks+=1
    print(dict(exact_core_fringe_checks=checks),flush=True)


def log_binom(n,k):
    if not 0<=k<=n: return float('-inf')
    return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/log(2)


def row_distribution(coverage):
    import numpy as np
    pmf=np.array([1.0])
    conditioning=0.0
    for count in coverage:
        p=count/M
        nxt=np.zeros(len(pmf)+1)
        nxt[:-1]+=pmf*(1-p);nxt[1:]+=pmf*p
        pmf=nxt
        if 0<count<M:
            conditioning-=log_binom(M,count)+count*log2(p)+(M-count)*log2(1-p)
    support=np.flatnonzero(pmf)
    return support,np.log2(pmf[support]),conditioning


def conditioned_rank(coverage,j,deaths,theta):
    """Positive conditional moment, minimized in total-coverage marker.

    Returns log2 bound for one mask before any support union. Actual
    regional quotas are retained. Projective normalization saves one
    field dimension; this diagnostic assumes a nonempty valid mask.
    """
    import numpy as np
    support,logs,cost=row_distribution(coverage)
    target=sum(coverage)/M
    tilt=32*theta*np.maximum(128-support,0)+logs
    def value(rho):
        terms=tilt+rho*support
        largest=float(terms.max())
        total=largest+log2(float(np.exp2(terms-largest).sum()))
        weights=np.exp2(terms-total)
        return total,float(weights@support)
    if len(support)==1:
        rho=0.0
    else:
        left=-128.0;right=128.0
        for _ in range(80):
            middle=(left+right)/2
            if value(middle)[1]<target:left=middle
            else:right=middle
        rho=(left+right)/2
    moment=cost-rho*sum(coverage)+M*value(rho)[0]
    bound=theta+32*theta*(j-deaths-K-1)+moment
    top=sum(sorted(coverage,reverse=True)[:128])
    direct=1+32*(j-deaths-top-1)
    return min(0.0,bound,direct),rho,cost


def earliest_birth_tail(coverage,region):
    """Pinned birth row: no M union. Exact Poisson-binomial in float."""
    import numpy as np
    pmf=np.array([1.0])
    for i,count in enumerate(coverage):
        if i==region: continue
        p=count/M;nxt=np.zeros(len(pmf)+1)
        nxt[:-1]+=pmf*(1-p);nxt[1:]+=pmf*p;pmf=nxt
    probability=float(pmf[128:].sum())
    return log2(probability) if probability else float('-inf')


def route_diagnostics():
    # Fixed geometry profiles. No average-quartet substitution is used:
    # these are conditional on the listed regional coverage counts.
    cases=[
        ('uniform_1_4',[M//4]*256, N//16, 0),
        ('uniform_1_2',[M//2]*256, N//6, 0),
        ('uniform_3_4',[3*M//4]*256, W, 0),
        ('uniform_15_16',[15*M//16]*256,W,0),
        ('six_update_bursts',[7*M//16]*256,7*N//16,N//16),
        ('seven_update_bursts',[8*M//16]*256,W,N//16),
        ('short_full_burst',[0]*140+[M]*116,W,512),
        ('sparse_prefix_full_suffix',[64]*140+[M]*116,W,0),
    ]
    for name,coverage,j,deaths in cases:
        choices=[(conditioned_rank(coverage,j,deaths,t)[0],t)
                 for t in (1/32,1/16,1/8,1/4,1/2,3/4,1)]
        bound,theta=min(choices)
        birth=next(i for i,v in enumerate(coverage) if v)
        nonzero_possible=sum(v>0 for v in coverage)>=129
        print(dict(profile=name,covered=sum(coverage),j=j,deaths=deaths,
                   minimum_nonzero_span=nonzero_possible,
                   log2_per_mask=round(bound,3),theta=theta,
                   log2_naive_support_union=round(bound+log_binom(N,j)+log_binom(j,deaths),3),
                   log2_pinned_birth=round(earliest_birth_tail(coverage,birth),3)),flush=True)


def uniform_route_moment(covered,theta):
    """Profile-uniform conditional census; total coverage stays fixed.

    Convex ordering replaces regional probabilities by their mean.
    The universal quota cost is used, not the fictitious equal-profile cost.
    """
    import numpy as np
    p=covered/N
    if p==1:return -32*theta*K
    assert 0<p<1
    support=np.arange(257,dtype=float)
    logs=np.array([log_binom(256,k) for k in range(257)])
    logs+=support*log2(p)+(256-support)*log2(1-p)
    tilt=logs+32*theta*np.maximum(128-support,0)
    target=covered/M
    def value(rho):
        terms=tilt+rho*support;largest=float(terms.max())
        total=largest+log2(float(np.exp2(terms-largest).sum()))
        return total,float(np.exp2(terms-total)@support)
    left=-128.;right=128.
    for _ in range(70):
        middle=(left+right)/2
        if value(middle)[1]<target:left=middle
        else:right=middle
    rho=(left+right)/2
    return min(-16*theta*covered,
               -32*theta*K+256*log2(M+1)-rho*covered+M*value(rho)[0])


def product_envelope(theta,logz,logx,logh,group_outputs=False):
    """All source weights and all destination supports; floating upper.

    Exact polynomial coefficients underlie this diagnostic. Rounding is
    not directed, so a negative number would still need a verified replay.
    """
    import numpy as np
    if not (-64<=logz<=0 and -64<=logx<=64 and -64<=logh<=64):return float('inf')
    def ladd(a,b):
        hi=max(a,b);return hi+log2(1+2**(-abs(a-b)))
    logm=32*theta+logz
    birth=ladd(0,logm+logh)
    survivor=(theta*ladd(0,32+logz/theta) if group_outputs else ladd(0,logm))
    active=ladd(logz,survivor+logh)
    local=np.array([4*ladd(0,logm+logh+logx)]+
                   [log2(comb(4,i))+4*logx+i*active+(4-i)*birth for i in range(1,5)])*log(2)
    coeff=np.array([0.])
    for _ in range(128):
        nxt=np.full(len(coeff)+4,-np.inf)
        for i in range(5):nxt[i:i+len(coeff)]=np.logaddexp(nxt[i:i+len(coeff)],coeff+local[i])
        coeff=nxt
    ratios=coeff/log(2)-np.array([log_binom(512,u) for u in range(513)])-logh*np.arange(513)
    index=int(np.argmax(ratios))
    return 4096*float(ratios[index])+max(0.,-512*logh),index


def global_diagnostics(group_outputs=False):
    """One coarse scalar potential, not a finite certificate or optimizer."""
    for fraction in (0.25,0.5,0.75,0.875,0.9375,1.0):
        covered=round(N*fraction);best=(float('inf'),None)
        for theta in (1/16,1/8,1/4,1/2,3/4,1):
            route=uniform_route_moment(covered,theta)
            def objective(v):
                logz,logx,logh=v
                got=product_envelope(theta,logz,logx,logh,group_outputs)
                if isinstance(got,float):return got
                return theta-32*theta+route-W*logz-covered*logx+got[0]
            point=[-32*theta-0.3,-8*theta,0.]
            value=objective(point)
            for step in (4.,2.,1.,0.5,0.25):
                for _ in range(3):
                    improved=False
                    for axis in range(3):
                        for sign in (-1,1):
                            candidate=point.copy();candidate[axis]+=sign*step
                            trial=objective(candidate)
                            if trial<value:point=candidate;value=trial;improved=True
                    if not improved:break
            if value<best[0]:best=(value,(theta,point,product_envelope(theta,*point,group_outputs)[1]))
        print(dict(global_coverage_fraction=fraction,all_history_log2_upper=round(best[0],3),
                   theta=best[1][0],markers=best[1][1],maximizing_incoming_weight=best[1][2],
                   outputs_grouped_before_moment=group_outputs),flush=True)


def balanced_expansion(theta,logz,x):
    """Balanced empty/full potential; exact quartets, floating evaluation."""
    import numpy as np
    a=2**(32*theta);z=2**logz;m=a*z;c=(1+2**(32+logz/theta))**theta
    linear=1-x*c
    radical=(linear*linear+4*m*x*x*z)**0.5
    h=(2*x*z/(radical+linear) if linear>=0 else (radical-linear)/(2*m*x))
    g=1+m*h*x;d=x*(1+m*h)/g
    assert 0<d<=1+1e-12
    local=np.log(np.array([1,4*d**3,6*d*d,4*d,1]))
    coeff=np.array([0.])
    for _ in range(128):
        nxt=np.full(len(coeff)+4,-np.inf)
        for j in range(5):nxt[j:j+len(coeff)]=np.logaddexp(nxt[j:j+len(coeff)],coeff+local[j])
        coeff=nxt
    values=[(log_binom(512,u)-coeff[u]/log(2))/min(u,512-u) for u in range(1,512)]
    kappa=max(0.,min(values));worst=values.index(min(values))+1
    return log2(g),log2(h),kappa,worst


def expansion_diagnostics():
    choices=[]
    for theta in (1/32,1/16,1/8,1/4,1/2,1):
        for offset in (-1.,-.3,0.,.3,1.):
            logz=min(0.,-32*theta+offset)
            for x in (.5,.625,.75,.8125,.875,.9375,.98):
                g,h,kappa,worst=balanced_expansion(theta,logz,x)
                choices.append((theta,logz,x,g,h,kappa,worst))
    for fraction in (.25,.5,.625,.75,.8):
        covered=round(fraction*N);route={t:uniform_route_moment(covered,t) for t in (1/32,1/16,1/8,1/4,1/2,1)}
        best=None
        for theta,logz,x,g,h,kappa,worst in choices:
            base=(theta-32*theta+route[theta]-W*logz-covered*log2(x)
                  +N*g+max(0.,-512*h)+256*kappa)
            needed=max(0.,(base+100)/kappa) if kappa else float('inf')
            if best is None or needed<best[0]:best=(needed,theta,logz,x,kappa,worst,base)
        print(dict(coverage_fraction=fraction,required_intermediate_mass_fraction=best[0]/N,
                   theta=best[1],logz=best[2],coverage_marker=best[3],kappa=best[4],
                   worst_input_weight=best[5],unrestricted_log2_upper=best[6],
                   certificate=False),flush=True)


def partition_kappa(bits):
    import numpy as np
    d=2**(-bits)
    local=np.log(np.array([1,4*d**3,6*d*d,4*d,1]))
    coeff=np.array([0.])
    for _ in range(128):
        after=np.full(len(coeff)+4,-np.inf)
        for j in range(5):after[j:j+len(coeff)]=np.logaddexp(after[j:j+len(coeff)],coeff+local[j])
        coeff=after
    return max(0.,min((log_binom(512,u)-coeff[u]/log(2))/min(u,512-u) for u in range(1,512)))


def joint_death_expansion_diagnostic():
    # Point profiles only. The expression covers all histories IN a profile,
    # but a list of points is not coverage of the full parameter domain.
    kappa={s:partition_kappa(s) for s in (1/8,1/4,1/2,1,2,4,8)}
    profiles=[(.1,0,.35),(.25,0,.625),(.25,.01,.625),(.25,.05,.625),
              (.4,0,.625),(.4,.01,.625),(.4,.05,.75),(.4,.1,.75),
              (.45,0,.75),(.45,.01,.75),(.45,.05,.8),(.5,.1,.8),
              (.6,.05,.8),(.7,.05,.8)]
    for af,df,ff in profiles:
        active=round(af*N);deaths=round(df*N);covered=round(ff*N)
        minority=min(active,N-active);excess=covered-active
        zero=max(0,active+deaths-W);survivors=active-deaths
        if zero>survivors:continue
        count=min(log_binom(N,active),log_binom(active+511,511)
                  +log_binom(active,deaths)+log_binom(N-active,deaths))
        partition=min(0.,min(s*(excess+512)-kap*max(0,minority-256) for s,kap in kappa.items()))
        best=min((count+t-32*t+32*t*(active-zero)+t*log_binom(survivors,zero)
                  +uniform_route_moment(covered,t)+partition,t)
                 for t in (1/32,1/16,1/8,1/4,1/2,3/4,1))
        print(dict(A=af,D=df,V=ff,minority=minority/N,log2_history_count=round(count,2),
                   log2_partition_tail=round(partition,2),log2_joint_upper=round(best[0],2),
                   theta=best[1],certificate=False),flush=True)


def expansion_rational_certificate():
    # P encodes a uniform quartet's live-cover expansion at marker1/4.
    local=(256,16,96,256,256);poly=[1]
    for _ in range(128):
        new=[0]*(len(poly)+4)
        for i,value in enumerate(poly):
            for j,coefficient in enumerate(local):new[i+j]+=value*coefficient
        poly=new
    failures=[]
    for u in range(513):
        if poly[u]*(1<<min(u,512-u))>(1<<1024)*comb(512,u):failures.append(u)
    assert not failures,failures
    # Independently check the free-terminal path count used by the proof.
    history_checks=0
    for steps in range(1,11):
        counts={}
        for word in range(1<<steps):
            previous=0;activity=0;deaths=0
            for t in range(steps):
                bit=(word>>t)&1
                activity+=bit;deaths+=previous*(1-bit);previous=bit
            key=(activity,deaths);counts[key]=counts.get(key,0)+1
        for activity in range(steps+1):
            for deaths in range(steps+1):
                expected=(comb(activity,deaths)*comb(steps-activity,deaths)
                          if deaths<=min(activity,steps-activity) else 0)
                assert counts.get((activity,deaths),0)==expected
                history_checks+=1
        for deaths in range((steps+1)//2+1):
            assert sum(v for (a,d),v in counts.items() if d==deaths)==(
                comb(steps+1,2*deaths+1) if 2*deaths+1<=steps+1 else 0)
    # Sum all state histories by the exact one-coordinate odd-binomial sum.
    delta=F(1,1<<15)
    # e^(2/3)<2 follows from the geometric Taylor bound, hence ln2>2/3.
    x=F(2,3)
    assert 1+x+x*x/(2*(1-x/3))<2
    # Z1=sum_r binom(4097,2r+1)*delta^(2r)
    # <=4097*exp(4096^2*delta^2/6), since (2r+1)!>=6^r*r!.
    assert 512*4096**2*delta**2/6==F(4,3)
    # 4097^512 <2^6144*exp(1/8); exp(35/24)<2^3.
    assert F(1,8)+F(4,3)==F(35,24)<3*F(2,3)
    assert -31+1024+6147==7140
    assert 7140-7216==-76
    print(dict(exact_quartet_expansion_comparisons=513,
               exact_free_terminal_history_comparisons=history_checks,
               forbidden_class='minority_mass+32*z-ceil_log2_binom(I,z) >=32*shortened_outer_dimension+7216',
               all_history_failure_below_2_to=-76,full_code_certificate=False),flush=True)


def footprint_expansion_certificate():
    """Reverse-hole expansion, with an output-assisted footprint class."""
    expansion_rational_certificate()
    checked=0;full=15
    for first in range(16):
        for second in range(16):
            for third in range(16):
                covers=(0,first,second,third)
                holes=tuple(full^v for v in covers)
                covered=sum(v.bit_count() for v in covers)
                allowed=third.bit_count()+sum(covers[t].bit_count()
                    for t in (1,2) if covers[t+1]==full)
                lost=covered-allowed;expansion=0;omitted=0;deaths=0
                for t in (1,2,3):
                    live=full if holes[t] else 0
                    omitted+=(live&~holes[t-1]).bit_count()
                    expansion+=live.bit_count()-holes[t].bit_count()
                    deaths+=(holes[t]&~holes[t-1]).bit_count()
                assert lost==omitted and lost>=deaths
                assert lost>=expansion-third.bit_count()
                assert deaths==sum((covers[t]&~covers[t+1]).bit_count() for t in (1,2))
                checked+=1
    q=1<<32
    assert F(35,24)+F(2*N,q)<3*F(2,3)
    assert F(35,24)+1+F(N,q)<4*F(2,3)
    assert -31+1024+6148-7217==-76
    print(dict(exact_reverse_hole_geometry_checks=checked,
        footprint_class='mV+11*max(0,L-W) >=32*dshort+7217',
        all_history_failure_below_2_to=-76,full_code_certificate=False),flush=True)


def remainder_profile_diagnostic():
    """Coarse scalar profiles only: exposes losses, never proves coverage."""
    kappa={s:partition_kappa(s) for s in (1/8,1/4,1/2,1,2,4,8)}
    thetas=(1/32,1/16,1/8,1/4,1/2,3/4,1)
    route_cache={};worst=[];checked=0
    for mf in (.0005,.002,.01,.025,.05,.1,.2,.4):
        minority=round(mf*N);best_failure=(-float('inf'),None)
        for af in (.1,.25,.4,.44,.45,.46,.48,.5,.55,.65,.75):
            active=round(af*N)
            if minority>min(active,N-active):continue
            for df in (0.,.0005,.002,.01,.025,.05,.1,.2):
                deaths=round(df*N)
                if deaths>min(active,N-active):continue
                zero=max(0,active+deaths-W);survivors=active-deaths
                if zero>survivors:continue
                count=min(log_binom(N,active),
                    min(log_binom(active+511,511),log_binom(N-active+511,511))
                    +log_binom(active,deaths)+log_binom(N-active,deaths),
                    4096+log_binom(N,minority))
                for extra in (0.,.25,.5,1.,2.,3.):
                    covered=active+deaths+round(extra*minority)
                    if not 0<covered<53*N/64:continue
                    partition=min(0.,min(s*(covered-active+512)-kap*max(0,minority-256)
                                          for s,kap in kappa.items()))
                    values=[]
                    for theta in thetas:
                        key=(covered,theta)
                        if key not in route_cache:route_cache[key]=uniform_route_moment(covered,theta)
                        values.append((count+theta-32*theta+32*theta*(active-zero)
                            +theta*log_binom(survivors,zero)+route_cache[key]+partition,theta))
                    bound,theta=min(values);checked+=1
                    if bound>best_failure[0]:best_failure=(bound,(af,df,extra,covered/N,theta,count,partition))
        worst.append(best_failure)
        print(dict(minority_fraction=mf,worst_sample_log2=best_failure[0],
                   profile=best_failure[1],profiles_checked=checked,certificate=False),flush=True)


def main():
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20 GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    if '--expansion-cert' in sys.argv:expansion_rational_certificate()
    elif '--footprint-expansion-cert' in sys.argv:footprint_expansion_certificate()
    elif '--remainder' in sys.argv:remainder_profile_diagnostic()
    elif '--joint' in sys.argv:joint_death_expansion_diagnostic()
    elif '--expansion' in sys.argv:expansion_diagnostics()
    elif '--global' in sys.argv:global_diagnostics('--group-outputs' in sys.argv)
    else:
        exact_transfer_checks()
        exact_transfer_checks(True)
        rational_coverage_certificate()
        exact_fringe_checks()
        route_diagnostics()
    print(dict(elapsed_seconds=round(time.monotonic()-start,3),certificate=False),flush=True)
    timer.cancel()


if __name__=='__main__':main()
