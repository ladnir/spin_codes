"""Exact all-row-count exclusion retaining A-C_T rather than dropping C_T.

No numerical search, files, or construction changes. This is a restricted
activity/coverage statement over J, not a minimum-distance certificate.
"""
import math
from fractions import Fraction as F
from itertools import product

N=1<<21
ROWS=8192
FIELD_BITS=32
SLACK=65794


def remaining_rank_score_floor(active_rows):
    assert 1<=active_rows<=ROWS
    return 128*active_rows-SLACK+1


def certificate():
    # The row-set and support unions are complete power-set sums. There is
    # no additional factor ROWS: it is already included in 2^ROWS.
    exponent=1+N+ROWS-FIELD_BITS*SLACK
    assert exponent==-63
    checks=0
    for rows in range(1,9):
        assert sum(math.comb(rows,n) for n in range(1,rows+1))==2**rows-1
        for positions in range(1,9):
            assert sum(math.comb(positions,a) for a in range(positions+1))==2**positions
            assert (2**rows-1)*2**positions<2**(rows+positions)
            checks+=1
    assert remaining_rank_score_floor(6092)==713983
    print(dict(score='total state activity minus covered inactive positions',
               all_row_counts=True,all_row_sets=True,arbitrary_histories=True,
               conditional_on_any_fixed_geometry=True,free_terminal=True,
               failure_less_than='2^-63',remaining_floor_at_active_6092=713983,
               exact_power_set_checks=checks,arithmetic='integer',
               whole_code_distance=False,files_written=0),flush=True)
    return SLACK


def markov_certificate():
    """Normalize the universal support prior by exact small path sums.

    General proof integrates two Bernoulli transition parameters shared by
    all coordinate chains. Their Beta(1,1) prior has total mass one.
    """
    checks=0
    for chains,length in ((1,4),(2,3),(3,2)):
        size=chains*length;total=F(0);fixed_total=F(0)
        for bits in product((0,1),repeat=size):
            a=sum(bits);deaths=0;terminal=0
            for j in range(chains):
                previous=0
                for t in range(length):
                    current=bits[j*length+t];deaths+=previous*(1-current);previous=current
                terminal+=previous
            births=deaths+terminal;n0=size-a+terminal;n1=a-terminal
            mu=F(1,(n0+1)*math.comb(n0,births)*(n1+1)*math.comb(n1,deaths))
            total+=mu
            b=F(1,9);d=F(1,6)
            fixed=(1-b)**(size-a-deaths)*b**births*d**deaths*(1-d)**(a-deaths-terminal)
            fixed_total+=fixed
            assert mu>=fixed/F((size+1)**2);checks+=1
        assert total==fixed_total==1
    assert (N+1)**2<2**44
    for rows in range(1,9):
        assert sum(F(math.comb(rows,n),(rows+1)*math.comb(rows,n)) for n in range(rows+1))==1
    print(dict(universal_markov_prior_normalizations=3,exact_path_comparisons=checks,
               universal_row_prior_normalizations=8,
               all_row_counts=True,all_histories=True,free_terminal=True,
               conditional_on_any_fixed_geometry=True,failure_at_most='2^-74',
               simultaneous_all_transition_parameters=True,
               whole_code_distance=False,arithmetic='rational',files_written=0),flush=True)


def markov_weights(active_rows,birth=F(1,9),death=F(1,6)):
    """Upper coefficients give a weaker necessary linear score >=0.

    The exact universal-prior exclusion holds simultaneously for every
    birth/death choice. All counts multiplying rounded coefficients are
    nonnegative. FLOAT conversion is only for an untrusted screen.
    """
    from quartet_two_update import log_upper
    assert 0<birth<1 and 0<death<1 and 1<=active_rows<=ROWS
    stay0=1-birth;stay1=1-death
    la=log_upper(1,stay0)+log_upper(-1,stay1)
    ld=log_upper(1,stay0)+log_upper(1,stay1)+log_upper(-1,birth)+log_upper(-1,death)
    lt=log_upper(1,stay1)+log_upper(-1,birth)
    row_cost=log_upper(1,ROWS+1)+log_upper(1,math.comb(ROWS,active_rows))
    constant=-128*active_rows+F(log_upper(-N,stay0)+row_cost+1+44+74,32)
    return dict(activity=1+la/32,covered=F(-1),deaths=ld/32,terminal=lt/32,constant=constant)


def selected_prior_certificate():
    """Count supports and selected intersection zeros in ONE normalized prior.

    A zero state has two next outcomes; a live state has three, including
    the selected/unselected distinction on a live-to-live transition.
    Independent Dirichlet(1,1) and Dirichlet(1,1,1) rows are shared by all
    coordinate chains. Tests enumerate all decorated paths, not samples.
    """
    checks=0;normalizations=0;linear_checks=0
    beta=F(1,10);death=F(1,10);selected=F(1,6)
    stay=1-death-selected
    for chains,length in ((1,4),(2,3),(3,2)):
        size=chains*length;total=F(0);fixed_total=F(0)
        for bits in product((0,1),repeat=size):
            activity=sum(bits);deaths=0;terminal=0
            for j in range(chains):
                previous=0
                for t in range(length):
                    current=bits[j*length+t]
                    deaths+=previous*(1-current);previous=current
                terminal+=previous
            n0=size-activity+terminal;n1=activity-terminal
            births=deaths+terminal;intersections=activity-deaths-terminal
            for zeros in range(intersections+1):
                unselected=intersections-zeros
                multinomial=math.comb(n1,deaths)*math.comb(n1-deaths,zeros)
                denominator=(n0+1)*math.comb(n0,births)*math.comb(n1+2,2)*multinomial
                mu=F(1,denominator)
                fixed=(1-beta)**(size-activity-deaths)*beta**births
                fixed*=death**deaths*stay**unselected*selected**zeros
                multiplicity=math.comb(intersections,zeros)
                total+=multiplicity*mu;fixed_total+=multiplicity*fixed
                overhead=(n0+1)*math.comb(n1+2,2)
                assert mu>=fixed/overhead;checks+=1
            for output_cap in range(size+1):
                zeros=max(0,activity+deaths-output_cap)
                if zeros>intersections:continue
                # If the selected-zero coefficient is nonpositive, using
                # A+D-W instead of max(0,A+D-W) gives an UPPER score.
                for coefficient in (F(-1),F(-1,2),F(0)):
                    assert coefficient*zeros<=coefficient*(activity+deaths-output_cap)
                    linear_checks+=1
        assert total==fixed_total==1;normalizations+=1
    assert (N+2)**3<2**64
    print(dict(joint_prior_normalizations=normalizations,
               exact_decorated_path_comparisons=checks,
               exact_zero_elimination_checks=linear_checks,
               all_row_counts=True,all_histories=True,all_selected_zero_sets=True,
               free_terminal=True,conditional_on_any_fixed_geometry=True,
               simultaneous_all_transition_parameters=True,failure_at_most='2^-74',
               whole_code_distance=False,arithmetic='rational',files_written=0),flush=True)


def selected_prior_weights(active_rows,birth=F(1,10),death=F(1,10),selected=F(1,6),output_cap=943718):
    """Necessary linear score for low-output histories, rounded UP.

    Eliminate the required selected-zero count z=max(0,A+D-W). Its score
    coefficient -1+log2(stay/selected)/32 must be nonpositive. Derive the
    exact combined coefficients BEFORE directed rounding, because W enters
    with the opposite sign. No assumptions on realized state transitions.
    """
    from quartet_two_update import log_upper
    assert 1<=active_rows<=ROWS and 0<birth<1
    assert death>0 and selected>0 and death+selected<1
    stay0=1-birth;stay1=1-death-selected
    assert stay1<=selected*(1<<32)
    row_cost=log_upper(1,ROWS+1)+log_upper(1,math.comb(ROWS,active_rows))
    activity=(log_upper(1,stay0)+log_upper(-1,selected))/32
    deaths=-1+(log_upper(1,stay0)+log_upper(2,stay1)
               +log_upper(-1,birth)+log_upper(-1,death)+log_upper(-1,selected))/32
    terminal=(log_upper(1,stay1)+log_upper(-1,birth))/32
    constant=output_cap-128*active_rows
    constant+=(log_upper(-N,stay0)+log_upper(-output_cap,stay1)
               +log_upper(output_cap,selected)+row_cost+1+64+74)/32
    return dict(activity=activity,covered=F(-1),deaths=deaths,terminal=terminal,constant=constant)


def conditional_weighted_cut_closure(slope=F(24,100),death_weight=F(74,100),offset=5000):
    """Exact implication at r=2100; the needed ensemble cut is NOT proved.

    On the already proved activity good event A>=743225. If geometry also
    obeys C+(74/100)D >= (24/100)A-B for every support history and row set
    of this size, the selected-prior good event yields a contradiction for
    the reported B. A single sampled geometry is not that ensemble event.
    """
    from quartet_zero_neighbors import certified_row_bounds
    rows=2100;active=ROWS-rows
    geometry=certified_row_bounds(rows)
    floor=geometry['activity_per_active_row']*active+1
    weights=selected_prior_weights(active)
    assert weights['activity']<slope and weights['deaths']+death_weight<=0
    upper=weights['constant']+(weights['activity']-slope)*floor
    upper+=max(F(0),weights['terminal'])*512
    budget=(-upper.numerator)//upper.denominator
    assert upper+offset<0
    print(dict(inactive_rows=rows,activity_floor=floor,
               assumed_cut=dict(slope=str(slope),death_weight=str(death_weight)),
               tested_B=offset,upper_remaining_score=str(upper+offset),
               conservative_maximum_integer_B=budget,
               contradiction_verified=True,arithmetic='directed rational',
               ensemble_cut_tail_proved=False,other_row_counts_proved=False,
               whole_code_distance=False,files_written=0),flush=True)
    return upper


def exclusion_budget():
    # Includes the optional stronger single-row-count activity checkpoint.
    exponents=(44,63,60,63,74,59,74,60)
    assert sum(1<<(74-bits) for bits in exponents)<1<<(74-43)
    assert (1<<(74-43))+(1<<(74-41))<1<<(74-40)
    print(dict(proved_exclusion_allowance_less_than='2^-43',
               sufficient_remaining_target='2^-41',remaining_target_proved=False,
               whole_code_distance=False,arithmetic='integer'),flush=True)


if __name__=='__main__':certificate()
