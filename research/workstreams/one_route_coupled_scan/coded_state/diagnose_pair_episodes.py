"""Episode-aware bounds for rejected-uniform SL2 pairs; not a whole-code cert.

Standard library only. Serial, source-only, <=60 seconds, small integer data.
Setup is fixed across messages. The zero-defect bound holds for every route.
"""
import argparse
import math
import time

M, b, N, FIELD_BITS = 8192, 512, 1 << 21, 32


def logceil(n):
    assert n > 0
    return (n-1).bit_length()


def volume_bits(length, weight, field_bits=FIELD_BITS):
    """Integer upper bound for log2 of a q-ary Hamming ball."""
    assert length>=0 and weight>=0 and field_bits>=1
    weight=min(length,weight);q=1<<field_bits
    if weight==length:return field_bits*length
    # binom(length,j)*q^j is increasing through weight under this condition.
    if q*(length-weight+1)<weight:return field_bits*length
    return min(field_bits*length,
               field_bits*weight+logceil(weight+1)+logceil(math.comb(length,weight)))


def zero_defect(R,V,permille,coarse=False):
    """No internal births, no one-sided incoming pairs, V active coordinates.

    Incoming states Q1..Q_(L-1) are counted; Q_L is free. The initial state
    before the R-region suffix is zero. Union over amplitudes and supports,
    NOT over messages chosen before setup.
    """
    assert 1<=R<=256 and 0<=V<=R*M-b and V%2==0
    L=R*M//b;W=permille*N//1000
    # Any proper support dies within nine distinct pairing phases. Thus at
    # most eight proper incoming states follow a prefix of full states.
    full_max=V//b
    prefix_batches=min(L,full_max+9)
    regions_free=(prefix_batches*b+M-1)//M
    regions_input=min(R,regions_free+1)  # possible uncounted terminal birth
    if V==0 or regions_input<=128:
        return dict(R=R,V=V,permille=permille,excluded=True,reason='RS support in at most128 regions')
    output_slots=min(R*M,(prefix_batches+1)*b)
    # F is supported in the prefix. The last input is fixed when its
    # incoming state is zero; if it is live it is already in this prefix.
    dimension=min(V//2,max(0,regions_free-128)*M)
    volume=min(FIELD_BITS*output_slots,output_slots+FIELD_BITS*W) if coarse else volume_bits(output_slots,W)
    bound=1+volume+FIELD_BITS*(dimension-V)
    return dict(R=R,V=V,permille=permille,excluded=False,
                max_full_states=full_max,prefix_batches=prefix_batches,
                output_slots=output_slots,free_regions=regions_free,
                outer_shortening_dimension=dimension,random_coordinates=V,
                log2_class_upper=bound,all_routes=True,whole_code_certificate=False)


def cover_zero_defects(permille=460):
    """All aligned suffixes, all even V: analytic majorant + interval replay."""
    began=time.monotonic();W=permille*N//1000;worst=-(1<<60);item=None;cells=0;excluded=0
    for R in range(129,257):
        # Fixed floor(V/b): prefix and dimension cap stay fixed. d(V)-V
        # is nonincreasing even when the V/2 cap is active. Only left is needed.
        for a in range(R*M//b):
            lo=max(2,a*b);hi=min((a+1)*b-2,R*M-b)
            if lo>hi:continue
            first=zero_defect(R,lo,permille,coarse=True)
            last=zero_defect(R,hi,permille,coarse=True)
            cells+=1
            if first['excluded']:
                assert last['excluded'];excluded+=1;continue
            assert not last['excluded'] and last['log2_class_upper']<=first['log2_class_upper']
            if first['log2_class_upper']>worst:worst=first['log2_class_upper'];item=(R,lo,hi)
        if time.monotonic()-began>30:raise TimeoutError('30s all-zero-defect cap')
    # This uniform bound is summed over output words directly, not over V.
    # There are at most 256 aligned suffix choices. Total density loss <2.
    analytic=9+N+FIELD_BITS*(W-N//2+24*b)
    assert worst+8<=analytic
    print(dict(kind='entire_zero_defect_family',permille=permille,
               count_intervals=cells,excluded_intervals=excluded,worst_interval=item,
               worst_class_log2_upper=worst,analytic_family_log2_upper=analytic,
               all_routes=True,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


def test():
    checked=0
    for bits in (1,2,3):
        q=1<<bits
        for n in range(1,25):
            for w in range(n+1):
                actual=sum(math.comb(n,j)*(q-1)**j for j in range(w+1))
                assert actual<=1<<volume_bits(n,w,bits)
                checked+=1
    # Exhaust the support geometry, independently of field-amplitude counting.
    fixtures=0
    for d in (1,2,3,4):
        width=1<<d;full=(1<<width)-1
        def closed(mask,phase):
            return all(((mask>>j)&1)==((mask>>(j^(1<<phase)))&1) for j in range(width))
        states=[[s for s in range(1<<width) if closed(s,p)] for p in range(d)]
        for start in states[0]:
            if start==full:continue
            reachable={start}
            for phase in range(1,d):
                reachable={v for v in states[phase] if any(v&~s==0 for s in reachable)}
            assert reachable=={0};fixtures+=1
    print(dict(exact_volume_checks=checked,exact_proper_support_starts=fixtures),flush=True)


def defect_test():
    """Exhaustive support-path DP, no sampling and no field assumptions."""
    d=2;width=1<<d;full=(1<<width)-1;length=20
    def singleton(mask,t):
        bit=1<<(t%d)
        return any(((mask>>j)&1)!=((mask>>(j^bit))&1) for j in range(width))
    states={(mask,0,int(0<mask<full)) for mask in range(full+1)}
    transitions=0
    for t in range(length-1):
        nxt=set()
        for mask,h,partial in states:
            for following in range(full+1):
                bad=int(singleton(mask,t) or bool(following&~mask))
                nxt.add((following,h+bad,partial+int(0<following<full)))
                transitions+=1
        states=nxt
    for mask,h,partial in states:
        h+=int(singleton(mask,length-1))
        assert partial<=d*h+d-1,(mask,h,partial)
    print(dict(exhaustive_support_paths_length=length,width=width,
               dynamic_states=len(states),transitions=transitions,
               verified='partial_batches <= d*h+d-1'),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--test',action='store_true')
    p.add_argument('--cover',action='store_true');p.add_argument('--defect-test',action='store_true')
    p.add_argument('--distance',type=int,default=460);args=p.parse_args()
    began=time.monotonic()
    if args.test:test()
    elif args.defect_test:defect_test()
    elif args.cover:cover_zero_defects(args.distance)
    else:
        for distance in (400,460,490,491,495):
            print(zero_defect(256,1048320,distance),flush=True)
            if time.monotonic()-began>60:raise TimeoutError('60s diagnostic cap; incomplete')
    print(dict(seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':main()
