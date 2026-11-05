"""Two-update, full-width flow envelope via a Gaussian moment majorant.

Construction unchanged. Bounded standard-library proof diagnostic. The
floating screen is not a distance certificate. No files or large arrays.
"""
import math
import shutil
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from random_pairing_demand_flow import causal_pair
from random_pairing_flow_kernel import kernel, falling
from random_pairing_flow_tail_screen import perron

N=512
M=N//2
P=F(1974,8192)
LAW=((1-P)**2,2*P*(1-P),P*P)


def matrix_product(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(5)) for j in range(5)] for i in range(5)]


def pair_data(z,y,v,law=LAW):
    pairs={};s=[[F(0) for _ in range(5)] for _ in range(5)]
    for a in range(5):
        for b in range(a,5):
            weights=[];marginal=[F(0)]*5
            for j,prob in enumerate(law):
                na,nb,_,loss=causal_pair(a-2,b-2,j);na+=2;nb+=2
                weight=prob*z**loss*y**j*v[na]*v[nb]/(v[a]*v[b])
                weights.append((na,nb,weight));marginal[na]+=weight/2;marginal[nb]+=weight/2
            total=sum(w for _,_,w in weights)
            s[a][b]=s[b][a]=total
            pairs[a,b]=(total,[x/total for x in marginal])
    return s,pairs


def positive_definite(a):
    """Exact LDL pivots; no floating eigenvalue acceptance."""
    residual=[row.copy() for row in a]
    for j in range(5):
        pivot=residual[j][j]
        if pivot<=0:return False
        for i in range(j+1,5):
            for k in range(j+1,5):
                residual[i][k]-=residual[i][j]*residual[j][k]/pivot
    return True


def majorize(s):
    # Entrywise larger and exactly positive definite. Coarse diagonal search
    # is deliberate: optimize the envelope only after this screen is useful.
    for numerator in range(101):
        extra=F(numerator,1000)
        candidate=[[s[i][j]+(extra if i==j else 0) for j in range(5)] for i in range(5)]
        if positive_definite(candidate):return candidate,extra
    raise AssertionError('No PSD majorant within the diagnostic search range')


def normalized_moment(c,q,m,exact=False):
    """E[(sum q_i G_i^2)^m]/(2m-1)!! for Cov(G)=c.

    Coefficient recurrence for det(I-t*diag(q)*c)^(-1/2), normalized by
    [t^m](1-t)^(-1/2). Exact arithmetic is used only for tiny validation.
    """
    cast=F if exact else float
    a=[[cast(q[i])*cast(c[i][j]) for j in range(5)] for i in range(5)]
    scale=max(sum(row) for row in a)
    a=[[x/scale for x in row] for row in a]
    power=[[cast(int(i==j)) for j in range(5)] for i in range(5)]
    traces=[]
    for _ in range(m):
        power=matrix_product(power,a);traces.append(sum(power[i][i] for i in range(5)))
    coeff=[cast(1)]
    denominator=cast(1)
    for k in range(1,m+1):
        coeff.append(sum(traces[i-1]*coeff[k-i] for i in range(1,k+1))/(2*k))
        denominator*=cast(2*k-1)/(2*k)
    return coeff[-1]/denominator*scale**m


def histogram_vector(histogram,v):
    return math.prod(v[a]**h for a,h in enumerate(histogram))


def exact_checks():
    # Four coordinates, all three tagged outcomes, actual uniform matching.
    v=tuple(F(x,100) for x in (100,94,89,86,85));z=F(27,25);y=F(7,6)
    law=(F(9,16),F(6,16),F(1,16))
    s,pairs=pair_data(z,y,v,law);c,extra=majorize(s)
    assert positive_definite(c)
    candidates={ab:total**2*normalized_moment(c,q,2,True) for ab,(total,q) in pairs.items()}
    bound=max(candidates.values())
    transition={};checked=0;worst=F(0)
    histograms=[]
    for a in range(5):
        for b in range(5-a):
            for d in range(5-a-b):
                for e in range(5-a-b-d):histograms.append((a,b,d,e,4-a-b-d-e))
    for h in histograms:
        transition[h]=kernel(h,False)
    one={h:sum(prob*z**loss*y**tags*histogram_vector(new,v)
               for (new,loss,tags),prob in transition[h].items()) for h in histograms}
    for h in histograms:
        direct=sum(prob*z**loss*y**tags*one[new] for (new,loss,tags),prob in transition[h].items())
        ratio=direct/histogram_vector(h,v)
        assert ratio<=bound;checked+=1;worst=max(worst,ratio)
    # Wick's fourth-moment identity independently checks the coefficient code.
    for _,q in pairs.values():
        aq=matrix_product([[q[i]*c[i][j] for j in range(5)] for i in range(5)],
                          [[q[i]*c[i][j] for j in range(5)] for i in range(5)])
        trace=sum(q[i]*c[i][i] for i in range(5))
        formula=(trace**2+2*sum(aq[i][i] for i in range(5)))/3
        assert normalized_moment(c,q,2,True)==formula
    print(dict(exact_all_histogram_two_step_checks=checked,exact_fourth_moment_checks=15,
               diagonal_majorant=float(extra),max_exact_ratio=float(worst),
               gaussian_bound=float(bound),not_a_distance_certificate=True),flush=True)


def homogeneous_envelope(z,y,deadline):
    """Necessary homogeneous-start envelope, not an all-histogram bound.

    Positive homogeneous operator iteration. Its min/max coordinate ratios
    bracket the optimal homogeneous-start envelope for fixed z,y. All
    values here are floating point, so the bracket is diagnostic only.
    """
    _,initial,_=perron(math.log(z),math.log(y))
    v=[math.sqrt(x) for x in initial]
    law=tuple(float(x) for x in LAW)
    for iteration in range(100):
        s,pairs=pair_data(z,y,v,law)
        ratios=[]
        for a in range(5):
            total,q=pairs[a,a]
            # Equal-state pairs output equal states. The formal Wick sum
            # is exact here even when s is not positive semidefinite.
            ratios.append((M*math.log(total)+math.log(normalized_moment(s,q,M)))/N)
        low,high=min(ratios),max(ratios)
        if high-low<1e-10:break
        v=[x*math.exp(value-high) for x,value in zip(v,ratios)]
        norm=max(v);v=[x/norm for x in v]
        assert time.monotonic()<deadline,'bounded necessary-envelope iteration'
    return low*N,high*N,v,iteration+1


def screen():
    z=F(27,25);y=F(7,6)
    _,eigenvector,_=perron(math.log(float(z)),math.log(float(y)))
    v=tuple(F(round(math.sqrt(x)*10**9),10**9) for x in eigenvector)
    s,pairs=pair_data(z,y,v);c,extra=majorize(s)
    cases=[];homogeneous=[];formal=[]
    for ab,(total,q) in pairs.items():
        moment=normalized_moment(c,q,M)
        cases.append((M*math.log(float(total))+math.log(moment),ab))
        unmajorized=normalized_moment(s,q,M)
        formal.append((M*math.log(float(total))+math.log(unmajorized),ab))
        if ab[0]==ab[1]:homogeneous.append(formal[-1])
    block_log,worst=max(cases)
    one_log=M*math.log(float(max(max(row) for row in s)))
    boundary=N*math.log(float(v[2]/min(v)))
    first=N*math.log(float(1-P+P*y))
    comb=math.lgamma(8193)-math.lgamma(1975)-math.lgamma(6219)
    regional=comb+1974*math.log(float(P))+6218*math.log(float(1-P))
    exponent=(2047*block_log+one_log+boundary+first-371901*math.log(float(z))
              -256*1974*math.log(float(y))+comb-256*regional)
    one_only=4095*one_log+boundary+first-371901*math.log(float(z))-256*1974*math.log(float(y))+comb-256*regional
    print(dict(width=N,diagonal_majorant=float(extra),worst_initial_pair=worst,
               two_step_log_per_pair=block_log/M,
               repeated_one_step_log_per_pair=2*one_log/M,
               full_selected_union_exponent_bits=exponent/math.log(2),
               one_step_selected_union_exponent_bits=one_only/math.log(2),
               exact_homogeneous_two_step_floor_bits=(exponent+2047*(max(homogeneous)[0]-block_log))/math.log(2),
               non_psd_formal_diagnostic_bits=(exponent+2047*(max(formal)[0]-block_log))/math.log(2),
               floating_screen_not_certificate=True),flush=True)
    tuned_floor(z,y,comb,regional)


def tuned_floor(z,y,comb,regional):
    low,high,tuned,iterations=homogeneous_envelope(float(z),float(y),time.monotonic()+18)
    # Necessary spectral diagnostic; nonnegative boundary costs omitted.
    centered_low=low-N*2*float(P)*math.log(float(y))
    centered_high=high-N*2*float(P)*math.log(float(y))
    selected_low=(2047*centered_low-371901*math.log(float(z))+comb-256*regional)/math.log(2)
    selected_high=(2047*centered_high-371901*math.log(float(z))+comb-256*regional)/math.log(2)
    print(dict(tuned_homogeneous_floor_bits=(selected_low,selected_high),
               endpoint_weights=tuned,iterations=iterations,
               fixed_markers_only=True,boundary_terms_omitted=True,
               floating_screen_not_certificate=True),flush=True)


def marker_ladder(deadline):
    comb=math.lgamma(8193)-math.lgamma(1975)-math.lgamma(6219)
    regional=comb+1974*math.log(float(P))+6218*math.log(float(1-P))
    best=(float('inf'),None)
    for theta in (.10,.14,.18):
        for ratio in (1.7,2.,2.3):
            eta=theta*ratio
            low,high,_,iterations=homogeneous_envelope(math.exp(theta),math.exp(eta),deadline)
            def exponent(value):
                return (2047*(value-2*N*float(P)*eta)-371901*theta+comb-256*regional)/math.log(2)
            candidate=(exponent(low),exponent(high))
            print(dict(theta=theta,eta=eta,homogeneous_floor_bits=candidate,iterations=iterations),flush=True)
            best=min(best,(candidate[1],(theta,eta)))
    print(dict(best_necessary_screen=best,all_histogram_certificate=False),flush=True)


def homogeneous_collision(a,n,m,r):
    """Exact E[sum h_next^2] from a homogeneous auxiliary state."""
    outputs=[]
    for j in range(3):
        left,right,_,_=causal_pair(a,a,j)
        assert left==right
        outputs.append(left)
    equal=F(0)
    for j in range(3):
        for k in range(3):
            if outputs[j]==outputs[k]:
                equal+=F(math.comb(2,j)*math.comb(2,k)*falling(r,j+k)*falling(m-r,4-j-k),falling(m,4))
    return 2*n+n*(n-2)*equal


def collision_feature():
    checked=0
    for a in range(-2,3):
        h=tuple(4*int(i==a+2) for i in range(5))
        law=kernel(h,True,(16,5))
        direct=sum(prob*sum(x*x for x in new) for (new,_,_),prob in law.items())
        assert direct==homogeneous_collision(a,4,16,5);checked+=1
    drift=[(a,float(F(1,100)*(homogeneous_collision(a,512,8192,1974)/512-512))) for a in range(-2,3)]
    print(dict(exact_regional_collision_checks=checked,log_potential_mean_changes=drift,
               potential='exp(sum h_a^2 / (100*n))',mean_only_not_exponential_bound=True),flush=True)


def main():
    assert shutil.disk_usage(Path.cwd()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    exact_checks()
    collision_feature()
    if '--ladder' in sys.argv:marker_ladder(began+22)
    elif '--test-only' not in sys.argv:screen()
    assert time.monotonic()-began<25
    print(dict(seconds=time.monotonic()-began))


if __name__=='__main__':main()
