"""Exact all-r cover for the allowed-input window-concentration event.

There are8192 source groups,4096 routing pairs per region, and256 regions.
The pair permutation in each region is independent and uniform. A window
of44 consecutive packets contains1408 pair slots. The event considered
here is: some set T of at least7 groups has ALL its allowed routing
positions inside one such window, separately in every region. The window
may differ between regions. This is a restriction on allowed positions,
not merely the nonzero entries of a particular outer word.

If the event holds for T, it holds for any seven-group subset. Fix one
such subset. It occupies at least4 routing pairs in either matching.
Select4 of those pairs deterministically before sampling the route.
Their destinations form a uniform four-subset of the4096 slots.

If four destinations lie in a contiguous1408-slot window, one destination
anchors a cyclic interval containing the other three among its next1407
slots. Allowing cyclic intervals only enlarges the event. A union over the
four anchors and sampling without replacement gives

    p4 <= 4 * C(1407,3)/C(4095,3)
        <= 4 * (1407/4095)^3 = 4 * (67/195)^3.

Independent regions and a union over seven-group subsets therefore give

    U = C(8192,7) * [4 * (67/195)^3]^256.

This covers every |T|>=7 simultaneously without a union over their sizes.
It does not cover arbitrary low-weight codewords: an outer word may be
zero on allowed positions outside a window. No encoder or prior proof
helper is modified. The inequality and maximal integer B with U<=2^-B
are verified using exact integer/Fraction arithmetic, not floating logs.

Run python -B. One stdout-only invocation;20GiB free-space guard,5-second
deadline,under16MiB estimated data. No files, caches, or numerical sweep.
"""

from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
DEADLINE=0.0


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("bounded5-second invocation exhausted")


def cover(group_count=8192,subset_size=7,regions=256,slots=4096,window_slots=1408):
    """Exact union bound; no use of a floating logarithm."""
    if not (group_count%2==0 and slots==group_count//2 and 1<=subset_size<=group_count
            and 1<=regions and 1<=window_slots<slots):
        raise ValueError("invalid paired-route geometry")
    anchors=(subset_size+1)//2
    one_region=min(Fraction(1),anchors*Fraction(window_slots-1,slots-1)**(anchors-1))
    union=min(Fraction(1),math.comb(group_count,subset_size)*one_region**regions)
    return anchors,one_region,union


def integer_margin(bound):
    if not 0<bound<=1:
        raise ValueError("requires a positive probability upper bound at most one")
    numerator,denominator=bound.numerator,bound.denominator
    bits=max(0,denominator.bit_length()-numerator.bit_length())
    while numerator<<bits>denominator:
        bits-=1
    while numerator<<(bits+1)<=denominator:
        bits+=1
    assert numerator<<bits<=denominator<numerator<<(bits+1)
    return bits


def self_test():
    anchor_cases=subsets=matching_cases=margin_cases=0
    for slots in (6,8,10):
        for active in (2,3,4):
            for width in range(active,slots):
                guard()
                successes=0
                for occupied in itertools.combinations(range(slots),active):
                    subsets+=1
                    has_window=any(all((x-start)%slots<width for x in occupied)
                                   for start in range(slots))
                    has_anchor=any(all((x-anchor)%slots<width for x in occupied)
                                   for anchor in occupied)
                    assert has_window==has_anchor
                    successes+=has_window
                probability=Fraction(successes,math.comb(slots,active))
                exact_anchor=active*Fraction(math.comb(width-1,active-1),math.comb(slots-1,active-1))
                ratio_anchor=active*Fraction(width-1,slots-1)**(active-1)
                assert probability<=exact_anchor<=ratio_anchor
                anchor_cases+=1
    for groups in (6,8):
        for size in range(3,groups+1):
            for occupied in itertools.combinations(range(groups),size):
                chosen=set(occupied[:3])
                for shift in (0,1):
                    pairs=[{(i+shift)%groups,(i+shift+1)%groups} for i in range(0,groups,2)]
                    all_pairs={i for i,pair in enumerate(pairs) if pair.intersection(occupied)}
                    sub_pairs={i for i,pair in enumerate(pairs) if pair.intersection(chosen)}
                    assert sub_pairs<=all_pairs and len(sub_pairs)>=2
                    matching_cases+=1
    for n in range(1,32):
        for d in range(n,64):
            value=Fraction(n,d)
            b=integer_margin(value)
            assert value<=Fraction(1,2**b) and value>Fraction(1,2**(b+1))
            margin_cases+=1
    return {"cyclic_anchor_parameter_cases":anchor_cases,"exhaustive_destination_subsets":subsets,
            "seven_subset_principle_toy_matching_checks":matching_cases,
            "exact_integer_margin_checks":margin_cases}


def main():
    global DEADLINE
    if len(sys.argv)!=1:
        raise ValueError("fixed certificate accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start=time.monotonic(); DEADLINE=start+5
    tests=self_test()
    anchors,regional,bound=cover()
    expected=math.comb(8192,7)*Fraction(4*67**3,195**3)**256
    assert anchors==4 and regional==Fraction(4*67**3,195**3) and bound==expected
    bits=integer_margin(bound)
    guard()
    print(json.dumps({"scope":"exact routing-event certificate only",
        "route_event_certificate":True,"whole_code_certificate":False,
        "covers_source_group_counts":"7 through8192 by seven-group subset reduction",
        "regions":256,"window_packets":44,"pair_slots_per_region":4096,
        "exact_bound":"C(8192,7) * (4 * (67/195)^3)^256",
        "maximal_integer_margin_bits":bits,
        "integer_comparison":{"numerator_times_2powB_le_denominator":True,
                              "numerator_times_2powBplus1_gt_denominator":True},
        "reduced_numerator_bits":bound.numerator.bit_length(),
        "reduced_denominator_bits":bound.denominator.bit_length(),
        "rational_sha256":hashlib.sha256(f"{bound.numerator}/{bound.denominator}".encode()).hexdigest(),
        "tests":tests,"seconds":time.monotonic()-start,"generated_bytes":0,
        "estimated_peak_bytes":16*(1<<20)},indent=2))


if __name__=="__main__":
    main()
