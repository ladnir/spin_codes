"""Joint RS/route/reset enumerator for the fixed three-shear research model.

Positive integer polynomial envelopes, no row independence assumption.
Serial standard-library checks, bounded caches, no generated files.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F
from functools import lru_cache

from diagnose_pair_defects import profile, ceildiv
from diagnose_pair_region_rank import early_joint, joint_cell

Q, M, L = 1 << 32, 8192, 4096


def multiply(a, b):
    out = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


@lru_cache(maxsize=512)
def episode_choices(length, aligned):
    """One closed episode: activation plus length incoming-state updates.

    Each regional fragment gives n+16z, with n=Q*ordinary_batches+length
    when it contains the single activation batch. Offsets are deterministic
    interval geometry, not another sampled permutation.
    """
    out = {}
    for offset in ((0,) if aligned else range(16)):
        remaining = length+1
        first = True
        poly = (1,)
        while remaining:
            take = min(16-offset, remaining)
            n = Q*(take-1)+length if first else Q*take
            poly = multiply(poly, (n,16))
            first = False
            remaining -= take
            offset = 0
        r = len(poly)-1
        previous = out.get(r, (0,)*len(poly))
        out[r] = tuple(max(x,y) for x,y in zip(previous,poly))
    return tuple(sorted(out.items()))


@lru_cache(maxsize=32)
def early_polynomials(h, a, c):
    began=time.monotonic()
    extra=h-c+a
    table={(0,0,0):(1,)}
    for index in range(c-1):
        nxt={}
        for (s,e,r),poly in table.items():
            for length in range(1,9+9*(extra-e)):
                cost=ceildiv(max(0,length-8),9)
                for count,other in episode_choices(length,index==0 and bool(a)):
                    key=s+length,e+cost,r+count
                    result=multiply(poly,other)
                    old=nxt.get(key,(0,)*len(result))
                    nxt[key]=tuple(max(x,y) for x,y in zip(old,result))
        table=nxt
        assert len(table)<20000
        if time.monotonic()-began>15:
            raise TimeoutError('15-second early-polynomial cap')
    grouped={}
    for (s,e,r),poly in table.items():
        grouped.setdefault((s,e),[]).append((r,poly))
    return grouped


@lru_cache(maxsize=4096)
def row_envelope(h,a,c,live,early,cost):
    """Per-row nonzero weighted enumerator, maximized over relaxed geometry."""
    tail=live-early+1  # terminal interval includes its activation, no reset
    if tail<=0:return F(0)
    complete,partial=divmod(tail,16)
    worst=F(0)
    for count,poly in early_polynomials(h,a,c).get((early,cost),()):
        if partial:
            poly=multiply(poly,(partial*Q,16))
            count+=1
        regions=complete+count
        assert regions<=256, 'split-region envelope needs <=256 nonzero slots'
        if regions<=128:continue
        maximum_zeros=regions-129
        numerator=0
        for j in range(maximum_zeros+1):
            for k in range(min(j,len(poly)-1)+1):
                removed_full=j-k
                if removed_full<=complete:
                    numerator += poly[k]*math.comb(complete,removed_full)*Q**(complete-removed_full)
        value=F(numerator,(16**count)*(Q**128))
        worst=max(worst,value)
    return worst


def reset_bound(h,a,c,live,early,cost):
    row=row_envelope(h,a,c,live,early,cost)
    if row==0:return -(1<<60)
    y=M*row
    if y>=1:return None
    # With no extra defects, all activations have pair-closed support. One
    # active RS row supplies at most one activation coordinate in a batch.
    minimum=2 if h-c+a==0 else 1
    value=(y*y/(1-y)) if minimum==2 else (2*y/(1-y))
    return value.numerator.bit_length()-value.denominator.bit_length()+1


def cover(h=6,seconds=50):
    began=time.monotonic()
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    bits=placement+9+80
    worst=-(1<<60); witness=None; base_cells=refined=weighted=0; residual=[]
    for a in (0,1):
        for c in range(1,h+a+1):
            for live in range(c,L):
                base=profile(h,a,c,live,460,True)
                base_cells+=1
                if base is None:continue
                bound=base['log2_upper']
                if h==6 and a==0 and c==6 and 1893<=live<=1908:
                    # Early proper intervals have length<=8. Their ten regions
                    # plus the terminal interval's <=119 regions must reach
                    # the RS minimum129, forcing the clean-crossing family.
                    bound=min(bound,crossing_per_placement(5))
                if bound+placement>-80:
                    local=-(1<<60);local_item=None
                    for early,cost,r in early_joint(h,a,c):
                        item=joint_cell(h,a,c,live,early,cost,r,bits)
                        refined+=1
                        if item is None:continue
                        value=bound if 'unresolved' in item else min(bound,item['upper'])
                        if h==6 and a==0 and c==6 and ceildiv(live-early+1,16)==119:
                            value=min(value,crossing_per_placement(5))
                        if value+placement>-80:
                            extra=reset_bound(h,a,c,live,early,cost)
                            weighted+=1
                            if extra is not None:value=min(value,extra)
                        if value>local:local,local_item=value,item
                        if time.monotonic()-began>seconds:
                            raise TimeoutError('50-second reset-weight cover cap; incomplete')
                    bound=local
                    if bound+placement>-80:residual.append(local_item)
                if bound>worst:worst,witness=bound,(a,c,live)
    print(dict(kind='no_full_reset_weighted_cover',defects=h,
               base_cells=base_cells,joint_cells=refined,weighted_cells=weighted,
               good_route_family_upper=worst+placement,route_failure_upper=-80,
               family_upper=max(worst+placement,-80)+1,worst=witness,
               residual_count=len(residual),residual_first=residual[:2],
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


@lru_cache(maxsize=64)
def crossing_per_placement(closed):
    factors={}
    for length in range(1,9):
        small=(length+1)//2
        large=length+1-small
        price=(32-(length-1).bit_length())//small
        factors[length]=F(small*large,256*(1<<price))
        for left in range(1,length+1):
            right=length+1-left
            actual_price=(32-(length-1).bit_length())//min(left,right)
            assert F(left*right,256*(1<<actual_price))<=factors[length]
    worst=-(1<<60)
    for maximum in range(1,9):
        per_row=factors[maximum]*max(factors[l] for l in range(1,maximum+1))**(closed-1)
        y=M*Q*per_row
        assert y<1
        probability=2*y**(1<<maximum)/(1-y)
        worst=max(worst,probability.numerator.bit_length()-probability.denominator.bit_length()+1)
    return worst


def clean_crossing_family():
    """5..64 clean proper closed episodes, then an unrestricted live tail.

    Each closed input interval crosses one regional boundary; touched
    regions of different intervals are disjoint. Total occupied regions=129.
    The final live interval ends at Q_(L-1), with no restriction on its
    internal births or one-sided pairs. This is not a whole-code family.
    """
    began=time.monotonic()
    factors={}
    for length in range(1,9):
        # Balanced splitting maximizes both route containment and the weaker
        # activation-rank bound. ceil(log2 length) loses <1 bit per death.
        small=(length+1)//2
        large=length+1-small
        price=(32-(length-1).bit_length())//small
        factors[length]=F(small*large,256*(1<<price))
        for left in range(1,length+1):
            right=length+1-left
            actual_price=(32-(length-1).bit_length())//min(left,right)
            assert F(left*right,256*(1<<actual_price))<=factors[length]
    worst=-(1<<60);witness=None;cases=0
    for closed in range(5,65):
        assert 129-2*closed>=1
        for maximum in range(1,9):
            per_row=factors[maximum]*max(factors[l] for l in range(1,maximum+1))**(closed-1)
            y=M*Q*per_row
            assert y<1
            min_rows=1<<maximum
            # At least min_rows eligible rows: backward matching closure
            # forces 2^length activation coordinates in each proper episode.
            # There are no unknown-coordinate choices in the forward bound.
            probability=2*y**min_rows/(1-y)
            per_place=probability.numerator.bit_length()-probability.denominator.bit_length()+1
            union=per_place+26*(closed+1)
            cases+=1
            if union>worst:worst,witness=union,(closed,maximum)
    final=worst+(cases-1).bit_length()
    assert max(crossing_per_placement(m)+26*(m+1) for m in range(5,65))==worst
    print(dict(kind='clean_crossing_reset_family',closed_episodes=(5,64),
               profile_cases=cases,worst=witness,worst_profile_upper=worst,
               entire_family_upper=final,terminal_internals='unrestricted',
               all_source_counts=True,route_included=True,
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)
    assert time.monotonic()-began<20


def test():
    began=time.monotonic()
    weighted=population=0
    # GF4 RS[4,2]: nonzero constants have no roots; each root is attained by
    # three nonconstant polynomials. Independently enumerate all coefficients.
    def mul(a,b):
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a&4:a^=7
        return out
    words=[tuple(u^mul(v,x) for x in range(4))
           for u,v in itertools.product(range(4),repeat=2) if u or v]
    for weights in itertools.product((F(0),F(1,4),F(1,2),F(3,4),F(1)),repeat=4):
        actual=sum((math.prod(weights[i] for i in range(4) if word[i]) for word in words),F(0))
        bound=16*math.prod(weights)+4*sum(math.prod(weights[i] for i in range(4) if i!=zero)
                                         for zero in range(4))
        assert actual<=bound;weighted+=1
        mean=sum(weights)/4
        for s in range(1,5):
            joint=sum(math.prod(weights[i] for i in subset)
                      for subset in itertools.combinations(range(4),s))/math.comb(4,s)
            assert joint<=mean**s;population+=1
    # Direct fragment products are bounded coefficientwise by the DP.
    dp=early_polynomials(2,0,2)
    geometries=0
    for length in range(1,9):
        for offset in range(16):
            batches=length+1;current=offset;poly=(1,);first=True
            while batches:
                take=min(16-current,batches)
                poly=multiply(poly,(Q*(take-1)+length if first else Q*take,16))
                first=False;batches-=take;current=0
            cap=dict(dp[length,0])[len(poly)-1]
            assert all(x<=y for x,y in zip(poly,cap));geometries+=1
    containment=closure=0
    for batches in range(1,5):
        for mask in range(1,1<<(4*batches)):
            coordinate_union=sum(any(mask>>(4*t+j)&1 for t in range(batches)) for j in range(4))
            assert mask.bit_count()<=batches*coordinate_union
            containment+=1
    for d in (1,2,3,4):
        width=1<<d;full=(1<<width)-1
        def closed(mask,phase):
            return all(((mask>>j)&1)==((mask>>(j^(1<<phase)))&1) for j in range(width))
        masks=[[s for s in range(1<<width) if closed(s,phase)] for phase in range(d)]
        for start in masks[0]:
            if start in (0,full):continue
            possible={start}
            for steps in range(1,d+1):
                if any(possible):assert start.bit_count()>=1<<steps
                closure+=1
                if steps<d:
                    possible={nxt for nxt in masks[steps] if any(nxt&~old==0 for old in possible)}
    print(dict(exact_weighted_rs_checks=weighted,weighted_without_replacement=population,
               fragment_envelope_checks=geometries,coordinate_containment_checks=containment,
               backward_closure_checks=closure,seconds=time.monotonic()-began),flush=True)
    assert time.monotonic()-began<20


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--test',action='store_true')
    p.add_argument('--cover',action='store_true')
    p.add_argument('--defects',type=int,default=6)
    p.add_argument('--clean-crossings',action='store_true')
    args=p.parse_args()
    assert 1<=args.defects<=8
    if args.test:test()
    elif args.clean_crossings:clean_crossing_family()
    elif args.cover:cover(args.defects)
    else:
        value=reset_bound(6,0,6,1893,5,0)
        print(dict(profile=(6,0,6,1893,5,0),per_placement_upper=value,
                   with_placement_union=value+186,whole_code_certificate=False),flush=True)


if __name__=='__main__':main()
