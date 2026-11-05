"""Exact quiet-constraint bound for consecutive256-group source supports.

Fix a consecutive set T of256 groups on the8192-cycle. Its two matching
profiles are (full,single)=(128,0),(127,2), in either order. Each profile
occurs in128 independently routed regions. Fix a44-packet window in each
region. Outside it, choose83 packets with both endpoint states zero,
omitting the final packet of every region. This also omits the sole global
packet without a feedback update. Any path whose endpoint-live packets
lie in the prescribed windows satisfies these selected quiet constraints.

For each selected packet, CX=0 has scalar rank sum_b min(4,ell_b).
Let D be the sum of these ranks across all selected packets. The established
exact-support scalar-shadow/outer-density interface bounds conditional
existence by min(1,eta*beta**256*q**(-D)), where q=256,
eta=q*q/(q-1)**2, beta=q**128*(q/(q-1))**32. This source statement concerns
messages with exact support T; it does not silently replace exact support
by containment. No low-output-weight condition is needed for this event.

Use y=3/4 and alpha=log(4/3)/log(256). Then 0<alpha<1/16 and
  (eta*beta**256)**alpha <= 2*y**(-32768)*(256/255)**512.
The region moment E[y**D] is an exact positive coefficient. One quiet
half-packet contains16 pair slots and contributes
  P(F,S)=sum_(f+s<=16) C(16,f)C(16-f,s) F**f S**s /2**s
          *sum_j C(s,j)y**(min(4,f+j)+min(4,f+s-j)).
There are166 quiet halves and1440 unconstrained pair slots per region.
The coefficient of F**f S**s in P**166*(1+F+S)**1440 is divided by
multinomial(4096;f,s,4096-f-s). Truncation to F<=128,S<=2 is exact because
all exponents are nonnegative. All retained P coefficients have the common
integer denominator2**18. The implementation uses positive integer
convolution only, and compares the final rational bound with powers of2.

The final union covers8192 consecutive source sets and85**256 choices of
non-wrapping44-packet windows. It covers this restricted path/support
class, not arbitrary supports or the whole code. No frozen helper is
imported or modified. One fixed stdout-only run,48-second deadline,
20GiB disk guard, estimated peak below64MiB; no files or caches.
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
DEN=1<<18


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("bounded48-second invocation exhausted")


def quiet_half(slots=16,rank=4,full_limit=128,single_limit=2):
    result=[[0]*(min(slots,full_limit)+1) for _ in range(single_limit+1)]
    for s in range(min(slots,single_limit)+1):
        for f in range(min(slots-s,full_limit)+1):
            value=0
            for j in range(s+1):
                exponent=min(rank,f+j)+min(rank,f+s-j)
                shift=18-s-2*exponent
                if shift<0:
                    raise ValueError("fixed denominator does not cover requested parameters")
                value+=math.comb(s,j)*3**exponent*(1<<shift)
            result[s][f]=math.comb(slots,f)*math.comb(slots-f,s)*value
    return result


def multiply(left,right,full_limit=128,single_limit=2):
    guard()
    width=min(full_limit+1,len(left[0])+len(right[0])-1)
    out=[[0]*width for _ in range(single_limit+1)]
    for ls,lrow in enumerate(left):
        for rs,rrow in enumerate(right[:single_limit-ls+1]):
            target=out[ls+rs]
            for i,a in enumerate(lrow[:width]):
                if a:
                    for j,b in enumerate(rrow[:width-i]):
                        if b:
                            target[i+j]+=a*b
                if i%16==0:
                    guard()
    return out


def positive_power(poly,exponent,full_limit=128,single_limit=2):
    result=[[1]]+[[0] for _ in range(single_limit)]
    while exponent:
        guard()
        if exponent&1:
            result=multiply(result,poly,full_limit,single_limit)
        exponent//=2
        if exponent:
            poly=multiply(poly,poly,full_limit,single_limit)
    return result


def extract_with_free(poly,full,single,free_slots):
    value=0
    for s in range(single+1):
        for f in range(min(full,len(poly[s])-1)+1):
            rf,rs=full-f,single-s
            if rf+rs<=free_slots:
                value+=poly[s][f]*math.comb(free_slots,rf)*math.comb(free_slots-rf,rs)
    return value


def margin(numerator,denominator):
    """Largest signed integer B satisfying numerator/denominator<=2**(-B)."""
    bits=denominator.bit_length()-numerator.bit_length()
    def le(b):
        return numerator<<b<=denominator if b>=0 else numerator<=denominator<<(-b)
    while not le(bits):
        bits-=1
    while le(bits+1):
        bits+=1
    assert le(bits) and not le(bits+1)
    return bits


def self_test():
    # Enumerate all half-slot labels: empty, full, left singleton, right
    # singleton. Fair singleton orientation has mass1/2, not a source count.
    checks=0
    for slots in (1,2,3):
        got=quiet_half(slots,1,slots,2)
        expected={}
        for labels in itertools.product(range(4),repeat=slots):
            f=labels.count(1); s=labels.count(2)+labels.count(3)
            if s>2:
                continue
            exponent=min(1,f+labels.count(2))+min(1,f+labels.count(3))
            expected[f,s]=expected.get((f,s),Fraction(0))+Fraction(3,4)**exponent/Fraction(2)**s
        for s,row in enumerate(got):
            for f,value in enumerate(row):
                assert Fraction(value,DEN)==expected.get((f,s),0)
                checks+=1
    # Explicit sequential products and dyadic exponentiation must agree,
    # including singleton degree2 and the unconstrained-slot convolution.
    poly=quiet_half(2,1,4,2)
    seq=[[1],[0],[0]]
    for exponent in range(1,6):
        seq=multiply(seq,poly,4,2)
        assert seq==positive_power(poly,exponent,4,2)
        for f in range(5):
            for s in range(3):
                direct=0
                for labels in itertools.product(range(3),repeat=2):
                    rf=labels.count(1); rs=labels.count(2)
                    if rf<=f and rs<=s and f-rf<len(seq[s-rs]):
                        direct+=seq[s-rs][f-rf]
                assert direct==extract_with_free(seq,f,s,2)
                checks+=1
    # The source-factor relaxation requires only these exact inequalities.
    assert Fraction(16,9)<2 and Fraction(256**2,255**2)<2
    for n in range(1,12):
        for d in range(1,12):
            b=margin(n,d)
            assert Fraction(n,d)<=Fraction(2)**(-b)
            assert Fraction(n,d)>Fraction(2)**(-b-1)
            checks+=1
    return checks


def main():
    global DEADLINE
    if len(sys.argv)!=1:
        raise ValueError("fixed exact check accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); DEADLINE=started+48
    checks=self_test()
    poly=positive_power(quiet_half(),166)
    numerators=[]; denominators=[]
    for f,s in ((128,0),(127,2)):
        numerators.append(extract_with_free(poly,f,s,1440))
        denominators.append(DEN**166*math.comb(4096,f)*math.comb(4096-f,s))
    guard()
    numerator=(numerators[0]*numerators[1])**128*2*4**32768*256**512*8192*85**256
    denominator=(denominators[0]*denominators[1])**128*3**32768*255**512
    bits=margin(numerator,denominator)
    guard()
    def integer_hash(value):
        return hashlib.sha256(value.to_bytes((value.bit_length()+7)//8,"big")).hexdigest()
    print(json.dumps({
        "scope":"exact consecutive256-support quiet-window class bound",
        "whole_code_certificate":False,"restricted_class_certificate":bits>0,
        "source_sets":8192,"regions":256,"window_packets":44,
        "window_choices_per_region":85,"selected_quiet_packets_per_region":83,
        "profiles":[[128,0],[127,2]],"regions_per_profile":128,
        "rank_mark":"3/4","source_factor_upper":"2*(4/3)^32768*(256/255)^512",
        "maximal_integer_margin_bits":bits,
        "integer_comparison":{"bound_le_2pow_minus_margin":True,"bound_gt_next_power":True},
        "regional_fraction_bits":[[n.bit_length(),d.bit_length()] for n,d in zip(numerators,denominators)],
        "unreduced_fraction_bits":[numerator.bit_length(),denominator.bit_length()],
        "numerator_sha256":integer_hash(numerator),"denominator_sha256":integer_hash(denominator),
        "exact_checks":checks,"elapsed_seconds":time.monotonic()-started,
        "generated_bytes":0,"estimated_peak_bytes":64*(1<<20),
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    },indent=2))


if __name__=="__main__":
    main()
