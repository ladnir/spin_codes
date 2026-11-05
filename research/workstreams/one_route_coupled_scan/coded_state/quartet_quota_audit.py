"""Directed obstruction to the Bernoulli moment route; exact quota reduction.

Proof-method obstruction only. Serial,90s watchdog,<32MiB integers,no files.
"""
from fractions import Fraction as F
from functools import lru_cache
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from quartet_two_update import log_interval,log_upper


@lru_cache(maxsize=24)
def choose_log(n,k):
    return log_interval(math.comb(n,k))


def linear_log(terms):
    lower=upper=F(0)
    for coefficient,(lo,hi) in terms:
        lower+=coefficient*(lo if coefficient>=0 else hi)
        upper+=coefficient*(hi if coefficient>=0 else lo)
    return lower,upper


def density_ratio(m,regions,r,unpaid):
    N=m*regions;R=r*regions;low,extra=divmod(unpaid,regions)
    assert 0<=unpaid<=R and low<=r
    terms=[(1,choose_log(N,unpaid)),(-1,choose_log(R,unpaid)),
           (-regions,choose_log(m,r)),(regions-extra,choose_log(m-low,r-low))]
    if extra:terms.append((extra,choose_log(m-low-1,r-low-1)))
    return linear_log(terms)


def density_tests():
    cases=0
    for m,regions,r in ((4,2,2),(4,3,2)):
        N=m*regions;R=r*regions
        for unpaid in range(R+1):
            low,extra=divmod(unpaid,regions)
            numerator=math.comb(m-low,r-low)**(regions-extra)
            if extra:numerator*=math.comb(m-low-1,r-low-1)**extra
            peak=F(numerator,math.comb(m,r)**regions*math.comb(R,unpaid))
            total=F(0);maximum=F(0)
            for subset in itertools.combinations(range(N),unpaid):
                counts=[sum(p//m==j for p in subset) for j in range(regions)]
                probability=F(math.prod(math.comb(m-u,r-u) if u<=r else 0 for u in counts),
                              math.comb(m,r)**regions*math.comb(R,unpaid))
                total+=probability;maximum=max(maximum,probability);assert probability<=peak;cases+=1
            assert total==1 and maximum==peak
            lo,hi=density_ratio(m,regions,r,unpaid)
            exact=peak*math.comb(N,unpaid)
            assert (lo,hi)==log_interval(exact) or lo<=log_interval(exact)[0]<=log_interval(exact)[1]<=hi
    print(dict(exhaustive_nested_quota_subsets=cases,normalization_and_balanced_maximum=True),flush=True)


def obstruction():
    N=1<<21;R=537600;c=111572;u=R-c;p=F(1523,2048)
    # For every nonnegative field multiplier, B/lambda<=Lmax. This uses
    # BOTH activity floors and c<=turnover_cap; it does not optimize numerically.
    Lmax=713983+c;a=F(39,512);constant=111127
    assert 743225<=Lmax and c-constant<=a*Lmax and c<=255463 and Lmax<N
    first=linear_log([(1,choose_log(R,c)),(Lmax,log_interval(p))])
    print(dict(covered=c,uniform_upper_on_B_over_lambda=Lmax,
               all_nonnegative_field_multipliers=True,all_positive_paid_markers=True,all_moment_bases_at_least_one=True,
               directed_log2_product_comparison_lower=str(first[0]),lower_float=float(first[0]),
               exact_spatial_temporal_moment_would_not_fix=True,whole_code_distance=False),flush=True)
    # Even retaining the unpaid count only by one Bernoulli evaluation and
    # dividing by its binomial point mass has a positive lower bound.
    second=linear_log([(1,choose_log(R,c)),(1,choose_log(N-Lmax,u)),(-1,choose_log(N,u))])
    print(dict(covered=c,directed_log2_single_binomial_conditioning_lower=str(second[0]),
               lower_float=float(second[0]),all_B_over_lambda_at_most=Lmax,
               all_Bernoulli_probabilities=True,exact_coefficient_extraction_not_ruled_out=True,
               whole_code_distance=False),flush=True)
    assert first[0]>0 and second[0]>0
    # At the previous frozen marker, even the ACTUAL unconditioned paid-set
    # polynomial already includes P=I; no Bernoulli replacement is needed.
    slope=F(707,2560);B=165556;logb=3
    frozen_lower=logb*(slope*N-B)-F(13,5)*(R-c)
    assert 6**5<2**13 and frozen_lower>0
    print(dict(frozen_marker='1/6',all_paid_endpoint_log2_lower=str(frozen_lower),
               lower_float=float(frozen_lower),bound_uses_log2_6_less_than_13_over_5=True,
               actual_paid_polynomial=True,not_a_fixed_count_coefficient=True),flush=True)


def exact_global_quota_obstruction():
    """A fixed-size suffix event rules out even the exact global coefficient.

    This event has zero probability under the nested regional quota law.
    No assertion about codewords or the actual construction follows.
    """
    N=1<<21;R=537600;c=111572;u=R-c;width=512
    Lmax=713983+c
    # Use a strict, complete-update suffix so H>=lambda*L>B for every
    # nonnegative multiplier, independently of the quartet partitions.
    L=width*(Lmax//width+1)
    assert L>Lmax and 0<=u<=N-L
    bound=linear_log([(1,choose_log(R,c)),(1,choose_log(N-L,u)),(-1,choose_log(N,u))])
    # Every complete region in the suffix must have all r positions paid.
    complete_regions=L//8192
    minimum_paid=complete_regions*2100
    assert minimum_paid>c
    assert bound[0]>38000
    print(dict(covered=c,suffix_positions=L,suffix_updates=L//width,
               uniform_fixed_count_event_probability='binom(N-L,u)/binom(N,u)',
               directed_log2_exact_global_quota_expression_lower=str(bound[0]),
               lower_float=float(bound[0]),all_moment_bases_at_least_one=True,
               all_current_nonnegative_field_multipliers=True,
               exact_coefficient_extraction_still_cannot_close=True,
               complete_regions_in_suffix=complete_regions,
               required_paid_positions_at_least=minimum_paid,
               actual_nested_quota_event_impossible=True,whole_code_distance=False),flush=True)
    # Small exact checks: count uniform subsets avoiding a suffix, and
    # count nested I,P pairs. The latter must vanish once full regions
    # in the suffix require more payments than allowed.
    cases=0
    for m,regions,r,c_small,L_small in ((4,3,2,1,4),(4,3,2,3,8),(4,2,1,0,4)):
        total=m*regions;unpaid=r*regions-c_small
        actual=sum(all(j<total-L_small for j in U)
                   for U in itertools.combinations(range(total),unpaid))
        assert actual==math.comb(total-L_small,unpaid)
        regional=[list(itertools.combinations(range(j*m,(j+1)*m),r)) for j in range(regions)]
        nested=0
        for blocks in itertools.product(*regional):
            I=tuple(itertools.chain.from_iterable(blocks))
            for P in itertools.combinations(I,c_small):
                cases+=1
                nested+=all(j<total-L_small or j in P for j in I)
        assert nested==0 and (L_small//m)*r>c_small
    print(dict(exact_nested_suffix_pairs_checked=cases,
               global_count_suffix_formula_checked=True,regional_impossibility_checked=True),flush=True)


def regional_free_boundary_obstruction():
    """The exact free-terminal regional product fails throughout onset three."""
    m=8192;regions=256;r=2100;width=512;c=111572
    floor_c,extra=divmod(c,regions)
    # In each region all unpaid positions lie before its last six updates.
    suffix=6*width;prefix=m-suffix;L=regions*suffix
    a=F(39,512);tau0=F(381,512);rank=713983;activity=743225
    cap=255463;constant=111127
    L0=rank+c;K0=a*L0-(c-constant)
    # B=lambda*L0-K0-nu*(rank+c-activity)-eta*(cap-c).
    # tau>=2lambda implies eta>=max(0,2lambda-tau0).
    # For lambda<=tau0/2, B/lambda<=L0-2K0/tau0.
    # Above it the same bound holds if tau0*(cap-c)-K0>=0.
    assert rank+c>=activity and K0>0 and tau0*(cap-c)-K0>=0
    ratio_max=L0-2*K0/tau0
    assert L>ratio_max
    terms=[(1,choose_log(regions,extra))]
    for count,multiplicity in ((floor_c,regions-extra),(floor_c+1,extra)):
        unpaid=r-count
        assert unpaid<=prefix
        terms.extend(((multiplicity,choose_log(r,count)),
                      (multiplicity,choose_log(prefix,unpaid)),
                      (-multiplicity,choose_log(m,unpaid))))
    lower,upper=linear_log(terms)
    assert lower>0
    slope=F(707,2560);tau=tau0;threshold=165556
    relaxed=slope*L
    # At an internal boundary the next quartet need not be active. Delete
    # the last selected update and pay the death of its predecessor.
    unrelaxed=relaxed-(slope+tau)*width*(regions-1)
    assert unrelaxed<threshold<relaxed
    print(dict(regional_suffix_updates=6,paid_low=floor_c,paid_high=floor_c+1,
               high_paid_regions=extra,directed_log2_free_boundary_product_lower=str(lower),
               lower_float=float(lower),all_current_multipliers_with_tau_at_least_twice_lambda=True,
               upper_B_over_lambda=str(ratio_max),chosen_total_suffix_positions=L,
               selected_block_score=str(relaxed),guaranteed_feasible_support_score=str(unrelaxed),
               selected_threshold=threshold,regional_quotas_exact=True,
               exact_paid_coefficient=True,no_read_degree_relaxation=True,
               no_scalar_paid_marker=True,whole_code_distance=False),flush=True)
    # Two-region blocks also fail at the selected field weights. Fix the
    # payment-count assignment and the long-suffix assignment, so nested
    # suffix events are not incorrectly counted as disjoint alternatives.
    long_regions=108
    terms=[]
    for count,updates,multiplicity in ((floor_c,6,regions-extra),
                                       (floor_c+1,6,extra-long_regions),
                                       (floor_c+1,7,long_regions)):
        unpaid=r-count
        terms.extend(((multiplicity,choose_log(r,count)),
                      (multiplicity,choose_log(m-width*updates,unpaid)),
                      (-multiplicity,choose_log(m,unpaid))))
    two_lower,_=linear_log(terms)
    two_L=L+long_regions*width
    two_score=slope*two_L-(slope+tau)*width*(regions//2)
    assert two_score>threshold and two_lower>0
    print(dict(regions_per_block=2,long_suffix_regions=long_regions,
               short_suffix_updates=6,long_suffix_updates=7,
               directed_log2_exact_two_region_product_lower=str(two_lower),
               lower_float=float(two_lower),guaranteed_block_score=str(two_score),
               selected_threshold=threshold,fixed_paid_and_suffix_assignments=True,
               exact_regional_quotas=True,exact_paid_coefficient=True,
               whole_code_distance=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        density_tests();obstruction();exact_global_quota_obstruction();regional_free_boundary_obstruction()
        for c in (27476,111572,255463):
            lower,upper=density_ratio(8192,256,2100,537600-c)
            print(dict(covered=c,unpaid=537600-c,directed_log2_density_lower=str(lower),
                       directed_log2_density_upper=str(upper),upper_float=float(upper),
                       exact_fixed_unpaid_count=True,uniform_subset_comparison=True,
                       regional_distribution_not_claimed_uniform=True,whole_code_distance=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
