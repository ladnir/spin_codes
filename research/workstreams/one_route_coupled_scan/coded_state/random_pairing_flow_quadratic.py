"""Finite-population mean and variational tests for a histogram potential.

Construction unchanged. Standard library, no files, bounded serial jobs.
Small exact tests validate the moment formula; screens are not certificates.
"""
import math
import shutil
import sys
import random
import time
from fractions import Fraction as F
from itertools import combinations_with_replacement, product, permutations
from pathlib import Path
from random_pairing_demand_flow import causal_pair
from random_pairing_flow_kernel import kernel, falling, matchings

TYPES=tuple(combinations_with_replacement(range(5),2))
QUADS=tuple(combinations_with_replacement(range(5),4))
P=F(1974,8192)
LAW=((1-P)**2,2*P*(1-P),P*P)
LOG_WEIGHT_NUMERATORS=(0,-84450914,-150416955,-196990597,-222300536)
# Separate proof candidate; the preserved certificates retain the tuple above.
RETUNED_LOG_WEIGHT_NUMERATORS=(0,-79200914,-154166955,-195740597,-227550536)
OUTPUT={}
LOSS={}
for a,b in TYPES:
    for j in range(3):
        na,nb,_,loss=causal_pair(a-2,b-2,j)
        OUTPUT[a,b,j]=tuple(int(c==na+2)+int(c==nb+2) for c in range(5))
        LOSS[a,b,j]=loss


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def counts(colors):return tuple(colors.count(c) for c in range(5))


def probability_colors(h,colors):
    multiplicities=counts(colors);degree=len(colors)
    numerator=math.factorial(degree)//math.prod(math.factorial(x) for x in multiplicities)
    numerator*=math.prod(falling(x,k) for x,k in zip(h,multiplicities))
    return numerator,falling(sum(h),degree)


def tables(pair_laws,joint_law=None):
    """Independent pair outcomes unless a regional joint tag law is supplied."""
    mean={};norm={};cross={}
    for a,b in TYPES:
        law=pair_laws[a,b]
        mean[a,b]=tuple(sum(law[j]*OUTPUT[a,b,j][c] for j in range(3)) for c in range(5))
        norm[a,b]=sum(law[j]*dot(OUTPUT[a,b,j],OUTPUT[a,b,j]) for j in range(3))
    for a,b,c,d in QUADS:
        total=0
        for ab,cd in (((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))):
            if joint_law is None:
                total+=dot(mean[ab],mean[cd])
            else:
                total+=sum(joint_law[j,k]*dot(OUTPUT[ab[0],ab[1],j],OUTPUT[cd[0],cd[1],k]) for j in range(3) for k in range(3))
        cross[a,b,c,d]=total/3
    return mean,norm,cross


def moments(h,prepared,exact=False):
    n=sum(h);assert n>=4 and n%2==0
    cast=F if exact else float
    mean,norm,cross=prepared
    next_mean=[cast(0)]*5;second=cast(0)
    for pair in TYPES:
        numerator,denominator=probability_colors(h,pair)
        chance=cast(numerator)/denominator
        for c in range(5):next_mean[c]+=(n//2)*chance*mean[pair][c]
        second+=(n//2)*chance*norm[pair]
    for quad in QUADS:
        numerator,denominator=probability_colors(h,quad)
        second+=(n//2)*(n//2-1)*cast(numerator)/denominator*cross[quad]
    return next_mean,second


def exact_tests():
    checked=0
    for regional in (False,True):
        if regional:
            m,r=16,5
            law=tuple(F(math.comb(2,j)*falling(r,j)*falling(m-r,2-j),falling(m,2)) for j in range(3))
            joint={(j,k):F(math.comb(2,j)*math.comb(2,k)*falling(r,j+k)*falling(m-r,4-j-k),falling(m,4))
                   for j in range(3) for k in range(3)}
        else:law=(F(9,16),F(6,16),F(1,16));joint=None
        prepared=tables({pair:law for pair in TYPES},joint)
        for colors in QUADS:
            h=counts(colors)
            next_mean,second=moments(h,prepared,True)
            actual=kernel(h,regional,(16,5) if regional else None)
            expected=[sum(prob*out[c] for (out,_,_),prob in actual.items()) for c in range(5)]
            expected_second=sum(prob*sum(x*x for x in out) for (out,_,_),prob in actual.items())
            assert next_mean==expected and second==expected_second
            checked+=1
    print(dict(exact_finite_population_moment_checks=checked),flush=True)


def grid(n=512,parts=8):
    assert n%parts==0
    for a in range(parts+1):
        for b in range(parts+1-a):
            for c in range(parts+1-a-b):
                for d in range(parts+1-a-b-c):
                    yield tuple(x*(n//parts) for x in (a,b,c,d,parts-a-b-c-d))


def screen():
    theta=.1;eta=.2;n=512
    # Frozen rational log weights, formerly selected by a spectral screen.
    u=[x/10**9 for x in LOG_WEIGHT_NUMERATORS]
    prepared=tables({pair:tuple(float(x) for x in LAW) for pair in TYPES})
    records=[]
    for h in grid():
        mean,second=moments(h,prepared)
        base=theta*float(LAW[0])*h[0]+dot(u,[x-y for x,y in zip(mean,h)])
        drift=second/n-sum(x*x for x in h)/n
        records.append((h,base,drift))
    comb=math.lgamma(8193)-math.lgamma(1975)-math.lgamma(6219)
    regional=comb+1974*math.log(float(P))+6218*math.log(float(1-P))
    allowance=(theta*371901-comb+256*regional)/4095
    print(dict(jensen_allowance_per_update=allowance,endpoint_costs_omitted=True),flush=True)
    for kappa in (0.,.005,.01,.02,.04):
        value,h=max((base+kappa*drift,h) for h,base,drift in records)
        print(dict(kappa=kappa,worst_grid_histogram=h,jensen_lower=value,
                   gap_to_allowance=value-allowance),flush=True)
    return records,u,theta,eta,allowance


def variational_reference(h,u,theta,eta,kappa):
    """A lower bound on the log moment via a changed pair-tag distribution.

    Pairings remain uniform. Fixed-point iteration only finds a proposal;
    the final variational inequality is valid for any positive proposal.
    """
    n=sum(h);law=tuple(float(x) for x in LAW)
    mean=[float(x) for x in h]
    for iteration in range(50):
        gradient=[u[c]+2*kappa*mean[c]/n for c in range(5)]
        proposals={}
        for a,b in TYPES:
            logits=[math.log(law[j])+theta*LOSS[a,b,j]+eta*j+dot(gradient,OUTPUT[a,b,j]) for j in range(3)]
            offset=max(logits);values=[math.exp(x-offset) for x in logits];mass=sum(values)
            proposals[a,b]=tuple(x/mass for x in values)
        updated,second=moments(h,tables(proposals))
        change=max(abs(x-y) for x,y in zip(updated,mean))
        mean=[(x+y)/2 for x,y in zip(updated,mean)]
        if change<1e-9:break
    mean,second=moments(h,tables(proposals))
    loss=tags=entropy=0.
    for pair in TYPES:
        numerator,denominator=probability_colors(h,pair)
        expected_pairs=(n//2)*numerator/denominator
        for j,chance in enumerate(proposals[pair]):
            loss+=expected_pairs*chance*LOSS[pair[0],pair[1],j]
            tags+=expected_pairs*chance*j
            entropy+=expected_pairs*chance*math.log(chance/law[j])
    value=theta*loss+eta*(tags-n*float(P))+dot(u,[x-y for x,y in zip(mean,h)])
    value+=kappa*(second-sum(x*x for x in h))/n-entropy
    return value,proposals,iteration+1


def variational(h,u,theta,eta,kappa):
    """Same proposal equation; compute second moments only after convergence.

    Each pair output has Euclidean diameter at most sqrt(8). Thus its
    tilted covariance is at most 2 I, and the mean map is 2*kappa-Lipschitz.
    Direct iteration is contractive for the frozen kappa=.08. The returned
    floating-point value is a screening lower bound, never a certificate.
    """
    assert 0<=kappa<.5
    n=sum(h);assert n>=4 and n%2==0
    law=tuple(float(x) for x in LAW)
    prepared=[]
    for a,b in TYPES:
        numerator,denominator=probability_colors(h,(a,b))
        expected=(n//2)*numerator/denominator
        outputs=tuple(OUTPUT[a,b,j] for j in range(3))
        logits=tuple(math.log(law[j])+theta*LOSS[a,b,j]+eta*j+dot(u,outputs[j]) for j in range(3))
        prepared.append(((a,b),expected,outputs,logits))
    mean=[float(x) for x in h]
    for iteration in range(50):
        proposals={};updated=[0.]*5
        for pair,expected,outputs,base in prepared:
            logits=[base[j]+2*kappa/n*dot(mean,outputs[j]) for j in range(3)]
            offset=max(logits);values=[math.exp(x-offset) for x in logits];mass=sum(values)
            proposal=tuple(x/mass for x in values);proposals[pair]=proposal
            for j,out in enumerate(outputs):
                weight=expected*proposal[j]
                for c in range(5):updated[c]+=weight*out[c]
        change=max(abs(x-y) for x,y in zip(updated,mean));mean=updated
        if change<1e-9:break
    mean,second=moments(h,tables(proposals))
    loss=tags=entropy=0.
    for (a,b),expected,_,_ in prepared:
        for j,chance in enumerate(proposals[a,b]):
            loss+=expected*chance*LOSS[a,b,j]
            tags+=expected*chance*j
            entropy+=expected*chance*math.log(chance/law[j])
    value=theta*loss+eta*(tags-n*float(P))+dot(u,[x-y for x,y in zip(mean,h)])
    value+=kappa*(second-sum(x*x for x in h))/n-entropy
    return value,proposals,iteration+1


def proposal_regression(deadline):
    """Compare proposals with the preserved damped solver, not just scores."""
    u=[x/10**9 for x in LOG_WEIGHT_NUMERATORS]
    points=list(grid(parts=4))+[(185,188,83,56,0),(103,103,102,102,102)]
    worst=0.;times=[0.,0.];iterations=[0,0]
    for h in points:
        results=[]
        for index,method in enumerate((variational_reference,variational)):
            start=time.monotonic();result=method(h,u,.1,.2,.08)
            times[index]+=time.monotonic()-start;iterations[index]+=result[2];results.append(result)
        old,new=results
        error=max(abs(old[1][pair][j]-new[1][pair][j]) for pair in TYPES for j in range(3))
        worst=max(worst,error)
        assert error<1e-10 and abs(old[0]-new[0])<1e-8,(h,old[0],new[0],error)
        assert time.monotonic()<deadline,'bounded proposal regression'
    print(dict(proposal_regression_points=len(points),max_probability_difference=worst,
               reference_seconds=times[0],direct_seconds=times[1],total_iterations=iterations),flush=True)


def variational_screen(records,u,theta,eta,allowance,deadline):
    for kappa in (.01,.02,.04,.08):
        selected=[h for _,h in sorted(((base+kappa*drift,h) for h,base,drift in records),reverse=True)[:8]]
        selected+= [tuple(512*int(i==a) for i in range(5)) for a in range(5)]
        best=(-float('inf'),None,0)
        for h in dict.fromkeys(selected):
            value,_,iterations=variational(h,u,theta,eta,kappa)
            best=max(best,(value,h,iterations))
            assert time.monotonic()<deadline,'bounded variational screen'
        print(dict(kappa=kappa,variational_lower=best[0],histogram=best[1],iterations=best[2],
                   gap_to_allowance=best[0]-allowance,selected_histograms_only=True),flush=True)


def deep_screen(records,u,theta,eta,allowance,deadline):
    kappa=.08;cache={}
    def evaluate(h):
        if h not in cache:cache[h]=variational(h,u,theta,eta,kappa)
        assert time.monotonic()<deadline,'bounded deep screen'
        return cache[h][0]
    best=max((evaluate(h),h) for h,_,_ in records)
    for step in (16,4,1):
        for _ in range(12):
            neighbors=[best]
            for a in range(5):
                if best[1][a]<step:continue
                for b in range(5):
                    if b==a:continue
                    h=list(best[1]);h[a]-=step;h[b]+=step;h=tuple(h)
                    neighbors.append((evaluate(h),h))
            updated=max(neighbors)
            if updated[0]<=best[0]+1e-12:break
            best=updated
    print(dict(kappa=kappa,grid_and_local_variational_lower=best[0],histogram=best[1],
               gap_to_allowance=best[0]-allowance,examined_histograms=len(cache),
               not_an_exhaustive_histogram_maximum=True),flush=True)
    return best[1],cache[best[1]][1],kappa


def importance_probe(h,proposals,u,theta,eta,kappa,deadline,samples=4096):
    n=sum(h);rng=random.Random(20261029)
    order=[a for a,count in enumerate(h) for _ in range(count)]
    lookup={}
    for a,b in TYPES:
        outcomes=[]
        for j in range(3):
            out=OUTPUT[a,b,j];indices=[c for c,x in enumerate(out) for _ in range(x)]
            score=theta*LOSS[a,b,j]+eta*j+dot(u,out)+math.log(float(LAW[j])/proposals[a,b][j])
            outcomes.append((indices[0],indices[1],score))
        lookup[a,b]=(proposals[a,b][0],sum(proposals[a,b][:2]),outcomes)
    constant=-eta*n*float(P)-dot(u,h)-kappa*sum(x*x for x in h)/n
    values=[]
    for sample in range(samples):
        rng.shuffle(order);out=[0]*5;score=constant
        for i in range(0,n,2):
            a,b=order[i],order[i+1]
            if a>b:a,b=b,a
            q0,q1,outcomes=lookup[a,b];uniform=rng.random()
            j=0 if uniform<q0 else (1 if uniform<q1 else 2)
            c,d,added=outcomes[j];out[c]+=1;out[d]+=1;score+=added
        score+=kappa*sum(x*x for x in out)/n;values.append(score)
        if sample%256==0:assert time.monotonic()<deadline,'bounded importance probe'
    center=sum(values)/samples
    ratios=[math.exp(x-center) for x in values];average=sum(ratios)/samples
    variance=sum((x-average)**2 for x in ratios)/(samples-1)
    print(dict(samples=samples,seed=20261029,importance_log_moment_estimate=center+math.log(average),
               estimated_log_standard_error=math.sqrt(variance/samples)/average,
               mean_log_weight=center,log_weight_range=(min(values),max(values)),
               numerical_diagnostic_not_upper_bound=True),flush=True)


def normal_quadrature(order):
    """Golub-Welsch via a small symmetric Jacobi eigensolver; diagnostic."""
    a=[[0.]*order for _ in range(order)];vectors=[[float(i==j) for j in range(order)] for i in range(order)]
    for i in range(order-1):a[i][i+1]=a[i+1][i]=math.sqrt(i+1)
    for _ in range(5000):
        _,p,q=max((abs(a[i][j]),i,j) for i in range(order) for j in range(i+1,order))
        if abs(a[p][q])<1e-14:break
        angle=.5*math.atan2(2*a[p][q],a[q][q]-a[p][p]);c=math.cos(angle);s=math.sin(angle)
        app,aqq,apq=a[p][p],a[q][q],a[p][q]
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq
        a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;a[p][q]=a[q][p]=0.
        for i in range(order):
            if i not in (p,q):
                oldp,oldq=a[i][p],a[i][q]
                a[i][p]=a[p][i]=c*oldp-s*oldq;a[i][q]=a[q][i]=s*oldp+c*oldq
            oldp,oldq=vectors[i][p],vectors[i][q]
            vectors[i][p]=c*oldp-s*oldq;vectors[i][q]=s*oldp+c*oldq
    else:raise AssertionError('quadrature eigenproblem')
    nodes=sorted((a[i][i],vectors[0][i]**2) for i in range(order))
    for power in range(2*order):
        expected=0 if power%2 else math.prod(range(1,power,2))
        actual=sum(w*x**power for x,w in nodes)
        magnitude=sum(w*abs(x)**power for x,w in nodes)
        assert abs(actual-expected)<1e-10*max(1,expected,magnitude)
    return nodes


def switch_range(f,active):
    result=0
    for ai,a in enumerate(active):
        for b in active[ai+1:]:
            differences=[f[x][a]-f[x][b] for x in active]
            result=max(result,max(differences)-min(differences))
    return result


def matching_variance_budget(f,h):
    """Range budget for exposing rare/extreme colors first.

    A phase choosing color a has at most h[a] exposures. Filling the
    earlier, larger range budgets first bounds every possible history.
    """
    active=tuple(a for a,count in enumerate(h) if count);cache={};best=None
    for priority in permutations(active):
        remaining=sum(h)//2-1;budget=0
        for index,a in enumerate(priority):
            if remaining<=0:break
            subset=tuple(sorted(priority[index:]))
            if subset not in cache:cache[subset]=switch_range(f,subset)**2
            steps=min(h[a],remaining);budget+=steps*cache[subset];remaining-=steps
        if best is None or budget<best:best=budget
    return best/8


def exact_range_tests():
    h=(2,2,2,2,0);colors=[a for a,count in enumerate(h) for _ in range(count)]
    f=[[F((a-b)**2+((a+b)%3),20) for b in range(5)] for a in range(5)]
    checked=0
    for priority in permutations(range(4)):
        ranks={a:i for i,a in enumerate(priority)};remaining=len(colors)//2-1;bound=F(0)
        for index,a in enumerate(priority):
            steps=min(h[a],remaining);bound+=steps*switch_range(f,priority[index:])**2/8;remaining-=steps
        for matching in matchings(tuple(range(len(colors)))):
            partner={}
            for a,b in matching:partner[a]=b;partner[b]=a
            live=set(range(len(colors)));cost=F(0)
            while len(live)>2:
                active=tuple(sorted(set(colors[x] for x in live)))
                cost+=switch_range(f,active)**2/8
                a=min(live,key=lambda x:(ranks[colors[x]],x));b=partner[a]
                live.remove(a);live.remove(b)
            assert cost<=bound;checked+=1
    print(dict(exact_range_budget_path_checks=checked),flush=True)


def verify_initial_rejection():
    """Exact Jensen rejection of the frozen kappa=.01 candidate."""
    h=(192,256,64,0,0);u=tuple(F(x,10**9) for x in LOG_WEIGHT_NUMERATORS)
    mean,second=moments(h,tables({pair:LAW for pair in TYPES}),True)
    lower=F(1,10)*LAW[0]*h[0]+dot(u,[x-y for x,y in zip(mean,h)])
    lower+=F(1,100)*(second-dot(h,h))/512
    choose=math.comb(8192,1974)
    assert choose>1<<6519
    probability=F(choose*987**1974*3109**6218,4096**8192)
    assert probability**2*(1<<13)<=1  # each region costs at least6.5 bits
    assert F(56,81)>F(69,100)  # positive atanh series lower bound for ln2
    exponent_lower=4095*lower-F(371901,10)+(6519+256*F(13,2))*F(69,100)
    assert exponent_lower>35
    print(dict(exact_initial_kappa_rejection=True,
               union_log_expression_exceeds_nats=35,code_counterexample=False),flush=True)


def matching_upper_quadrature(h,u,theta,eta,kappa,order,deadline):
    """Quadrature of a proved pointwise upper integrand, NOT certified integration."""
    n=sum(h);_,proposal,_=variational(h,u,theta,eta,kappa)
    center,_=moments(h,tables(proposal));beta=math.sqrt(2*kappa/n)
    base=[u[c]+2*kappa*center[c]/n for c in range(5)]
    constant=-eta*n*float(P)-dot(u,h)-kappa*(dot(h,h)+dot(center,center))/n
    pair_means={pair:(n//2)*probability_colors(h,pair)[0]/probability_colors(h,pair)[1] for pair in TYPES}
    # Four orthonormal contrasts; the common Gaussian direction cancels.
    basis=[[0.]*4 for _ in range(5)]
    for j in range(4):
        for i in range(j+1):basis[i][j]=1/math.sqrt((j+1)*(j+2))
        basis[j+1][j]=-(j+1)/math.sqrt((j+1)*(j+2))
    nodes=normal_quadrature(order);terms=[];active=[a for a,count in enumerate(h) if count]
    for index,point in enumerate(product(nodes,repeat=4)):
        coordinates=[x for x,_ in point];g=[beta*dot(row,coordinates) for row in basis]
        gradient=[x+y for x,y in zip(base,g)]
        f=[[0.]*5 for _ in range(5)];mean=0.
        for a,b in TYPES:
            logits=[math.log(float(LAW[j]))+theta*LOSS[a,b,j]+eta*j+dot(gradient,OUTPUT[a,b,j]) for j in range(3)]
            offset=max(logits);value=offset+math.log(sum(math.exp(x-offset) for x in logits))
            f[a][b]=f[b][a]=value;mean+=pair_means[a,b]*value
        budget=matching_variance_budget(f,h)
        value=constant+mean-dot(g,center)+budget
        terms.append(value+sum(math.log(w) for _,w in point))
        if index%256==0:assert time.monotonic()<deadline,'bounded Gaussian quadrature'
    maximum=max(terms)
    return maximum+math.log(sum(math.exp(x-maximum) for x in terms))


def endpoint_allowance(u,theta,eta,kappa,verbose=True):
    n=512;sorted_u=sorted(u)
    for size in range(1,6):
        threshold=(2*kappa+sum(sorted_u[:size]))/size
        if size==5 or threshold<=sorted_u[size]:break
    fractions=[max(0.,(threshold-x)/(2*kappa)) for x in u]
    assert abs(sum(fractions)-1)<1e-12
    minimum=dot(u,fractions)+kappa*dot(fractions,fractions)
    boundary=n*(u[2]+kappa-minimum)
    initial=n*math.log((1-float(P))*math.exp(-eta*float(P))+float(P)*math.exp(eta*(1-float(P))))
    comb=math.lgamma(8193)-math.lgamma(1975)-math.lgamma(6219)
    regional=comb+1974*math.log(float(P))+6218*math.log(float(1-P))
    allowance=(theta*371901-comb+256*regional-boundary-initial)/4095
    if verbose:
        print(dict(kappa=kappa,continuous_simplex_endpoint_cost=boundary,
                   ignored_initial_batch_cost=initial,allowance_with_endpoint_costs=allowance),flush=True)
    return allowance


def retune_potential(deadline):
    """Bounded minimax proposal search; lower-bound screens are NOT certificates.

    Only linear potential weights change. Validate each search stage on the
    full coarse grid and its neighboring maxima before the next stage.
    """
    theta=.1;eta=.2;kappa=.08
    original=tuple(x/10**9 for x in LOG_WEIGHT_NUMERATORS);u=original
    all_points=list(grid());points=set(all_points)
    points.add((185,188,83,56,0))
    def records(weights,candidates):
        allowance=endpoint_allowance(weights,theta,eta,kappa,False);result=[]
        for h in sorted(candidates):
            assert time.monotonic()<deadline,'bounded potential search'
            value,_,_=variational(h,weights,theta,eta,kappa)
            result.append((value-allowance,h))
        return sorted(result,reverse=True)
    ranked=records(u,points);active={h for _,h in ranked[:32]}
    print(dict(stage='original_lower_screen',worst=ranked[0],points=len(points)),flush=True)
    for step in (.001,.0005,.00025):
        for sweep in range(3):
            old=u;best=records(u,active)[0][0]
            for axis in range(1,5):
                for sign in (-1,1):
                    candidate=list(u);candidate[axis]+=sign*step;candidate=tuple(candidate)
                    score=records(candidate,active)[0][0]
                    if score<best:best=score;u=candidate
            if u==old:break
        ranked=records(u,points)
        for _,h in ranked[:8]:
            for size in (4,16):
                for a in range(5):
                    if h[a]<size:continue
                    for b in range(5):
                        if a==b:continue
                        neighbor=list(h);neighbor[a]-=size;neighbor[b]+=size
                        points.add(tuple(neighbor))
        ranked=records(u,points);active.update(h for _,h in ranked[:32])
        print(dict(stage='retuned_lower_screen',step=step,weights=u,worst=ranked[0],
                   points=len(points),active=len(active),not_an_upper_certificate=True),flush=True)
    # Keep both original and final values on precisely the same sample set.
    old=records(original,points);new=records(u,points)
    print(dict(original_worst_on_final_points=old[0],retuned_worst_on_final_points=new[0],
               final_weights=u,points=len(points),not_exhaustive=True),flush=True)
    for _,h in new[:3]:
        value=matching_upper_quadrature(h,u,theta,eta,kappa,5,deadline)
        print(dict(histogram=h,order5_upper_integrand_quadrature=value,
                   allowance=endpoint_allowance(u,theta,eta,kappa,False),not_directed=True),flush=True)


def gaussian_tail_screen(h,u,theta,eta,kappa):
    """Analytical radial tail formula; printed evaluation is floating point."""
    n=sum(h);_,proposal,_=variational(h,u,theta,eta,kappa)
    center,_=moments(h,tables(proposal));beta=math.sqrt(2*kappa/n)
    base=[u[c]+2*kappa*center[c]/n for c in range(5)]
    f=[[0.]*5 for _ in range(5)];mean=0.;gradient=[0.]*5
    for a,b in TYPES:
        logits=[math.log(float(LAW[j]))+theta*LOSS[a,b,j]+eta*j+dot(base,OUTPUT[a,b,j]) for j in range(3)]
        offset=max(logits);values=[math.exp(x-offset) for x in logits];mass=sum(values)
        f[a][b]=f[b][a]=offset+math.log(mass)
        numerator,denominator=probability_colors(h,(a,b));pairs=(n//2)*numerator/denominator
        mean+=pairs*f[a][b]
        for c in range(5):gradient[c]+=pairs*sum(values[j]/mass*OUTPUT[a,b,j][c] for j in range(3))
    constant=-eta*n*float(P)-dot(u,h)-kappa*(dot(h,h)+dot(center,center))/n+mean
    delta=switch_range(f,tuple(a for a,count in enumerate(h) if count))
    lipschitz=4*math.sqrt(2)*beta;steps=n//2-1
    # Hessian of each log gate moment <=2*beta^2 I; sum gives <=2*kappa I.
    residual=beta*math.sqrt(sum((x-y)**2 for x,y in zip(gradient,center)))
    a=.5-kappa-steps*lipschitz*lipschitz/8
    b=residual+steps*delta*lipschitz/4
    d=steps*delta*delta/8
    assert a>0
    for radius in (16,24):
        effective=a-b/radius;assert effective>0
        log_tail=constant+d+math.log1p(effective*radius*radius)-math.log(4*effective*effective)-effective*radius*radius
        print(dict(radius=radius,analytical_log_absolute_tail_bound=log_tail,
                   radial_quadratic_coefficient=a,numerical_evaluation_not_directed=True),flush=True)


def upper_grid(records,u,theta,eta,deadline):
    kappa=.08;cache={}
    def evaluate(h):
        if h not in cache:cache[h]=matching_upper_quadrature(h,u,theta,eta,kappa,3,deadline)
        return cache[h]
    best=max((evaluate(h),h) for h,_,_ in records)
    for step in (16,4,1):
        for _ in range(12):
            candidates=[best]
            for a in range(5):
                if best[1][a]<step:continue
                for b in range(5):
                    if b==a:continue
                    h=list(best[1]);h[a]-=step;h[b]+=step;h=tuple(h)
                    candidates.append((evaluate(h),h))
            updated=max(candidates)
            if updated[0]<=best[0]+1e-12:break
            best=updated
    allowance=endpoint_allowance(u,theta,eta,kappa)
    print(dict(upper_grid_examined=len(cache),order3_worst_histogram=best[1],
               order3_upper_integrand_quadrature=best[0],not_exhaustive_or_certified=True),flush=True)
    for order in (5,7,9):
        value=matching_upper_quadrature(best[1],u,theta,eta,kappa,order,deadline)
        print(dict(order=order,upper_integrand_quadrature=value,allowance=allowance,
                   selected_union_exponent_if_uniform=(value-allowance)*4095/math.log(2),
                   not_a_uniform_bound=True),flush=True)
    gaussian_tail_screen(best[1],u,theta,eta,kappa)


def retuned_upper_search(deadline):
    """Numerical upper-integrand search over the interior and every face.

    Quadrature is not a certificate. Preserve the fixed candidate weights;
    seek displaced maxima before spending time on directed coverage.
    """
    u=tuple(x/10**9 for x in RETUNED_LOG_WEIGHT_NUMERATORS)
    cache={}
    def evaluate(h):
        assert time.monotonic()<deadline,'bounded retuned upper search'
        if h not in cache:cache[h]=matching_upper_quadrature(h,u,.1,.2,.08,3,deadline)
        return cache[h]
    points=list(grid())+[(185,188,83,56,0),(189,168,115,40,0)]
    for h in points:evaluate(h)
    winners=[]
    for face in (None,0,1,2,3,4):
        best=max((value,h) for h,value in cache.items() if face is None or h[face]==0)
        for step in (16,4,1):
            for _ in range(8):
                choices=[best]
                for a in range(5):
                    if a==face or best[1][a]<step:continue
                    for b in range(5):
                        if b==face or a==b:continue
                        h=list(best[1]);h[a]-=step;h[b]+=step;h=tuple(h)
                        choices.append((evaluate(h),h))
                updated=max(choices)
                if updated[0]<=best[0]+1e-12:break
                best=updated
        winners.append((face,best))
        print(dict(face=face,order3_maximum=best[0],histogram=best[1],examined=len(cache)),flush=True)
    allowance=endpoint_allowance(u,.1,.2,.08,False)
    for h in dict.fromkeys(best[1] for _,best in winners):
        values=[matching_upper_quadrature(h,u,.1,.2,.08,order,deadline) for order in (5,7)]
        print(dict(histogram=h,orders=(5,7),numerical_upper_integrals=values,allowance=allowance,
                   not_exhaustive_or_directed=True),flush=True)
    print(dict(examined_histograms=len(cache),proof_weights=RETUNED_LOG_WEIGHT_NUMERATORS,
               whole_code_certificate=False),flush=True)


def main():
    assert shutil.disk_usage(Path.cwd()).free>=20*(1<<30),'disk guard'
    if '--retuned-upper-search' in sys.argv:
        began=time.monotonic();retuned_upper_search(began+85)
        print(dict(seconds=time.monotonic()-began,whole_code_certificate=False));return
    if '--proposal-regression' in sys.argv:
        began=time.monotonic();exact_tests();proposal_regression(began+24)
        print(dict(seconds=time.monotonic()-began,whole_code_certificate=False));return
    if '--retune-potential' in sys.argv:
        began=time.monotonic();retune_potential(began+50)
        print(dict(seconds=time.monotonic()-began,whole_code_certificate=False));return
    began=time.monotonic();exact_tests();exact_range_tests();verify_initial_rejection();args=screen()
    if '--upper-grid' in sys.argv:
        upper_grid(args[0],args[1],args[2],args[3],began+23)
    elif '--upper' in sys.argv:
        h=(185,188,83,56,0)
        allowance=endpoint_allowance(args[1],args[2],args[3],.08)
        for order in (5,7,9):
            value=matching_upper_quadrature(h,args[1],args[2],args[3],.08,order,began+23)
            print(dict(histogram=h,order=order,upper_integrand_quadrature=value,
                       gap_to_allowance=value-allowance,not_certified_quadrature=True),flush=True)
    elif '--deep' in sys.argv:
        h,proposal,kappa=deep_screen(*args,began+20)
        importance_probe(h,proposal,args[1],args[2],args[3],kappa,began+24)
    else:variational_screen(*args,began+22)
    assert time.monotonic()-began<25
    print(dict(seconds=time.monotonic()-began,no_exponential_upper_bound=True))


if __name__=='__main__':main()
