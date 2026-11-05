"""Exact consecutive-source size band for the44-packet quiet-window event.

This wrapper reuses the hash-pinned positive integer coefficient engine
from cheap_lane_quiet_window_cover.py. It checks every integer source size
r=7,...,256. An even size r=2m has matching profiles (m,0),(m-1,2);
an odd size r=2m+1 has profile (m,1) in both matchings. Each matching occurs
in128 regions. Every size has8192 distinct consecutive supports.

The path restriction is unchanged: in each region, all packets having a
nonzero endpoint state lie in some non-wrapping44-packet window. Choosing
83 quiet packets outside that window, excluding the last packet of every
region, gives the same rank moment table for every size. There are85
possible windows per region. This is a restricted path/support class,
not a distance certificate for arbitrary supports or the whole code.

Write q=256,y=3/4,alpha=log(4/3)/log(q),beta=q**128*(q/(q-1))**32.
The exact-support shadow prefactor has eta=q*q/((q-1)**2*(q-r+1)).
For r<=256, eta<=256**2/255**2<2 and 0<alpha<1/16, hence
  (eta*beta**r)**alpha <=2*y**(-128*r)*(256/255)**(2*r).
All remaining computations use exact integers. Every per-size bound is
rounded upward to a power of2. Their tail sums use one common dyadic
denominator, avoiding a large rational least-common multiple. The reported
tail is the smallest tail certified at40 bits by these dyadic ceilings;
no monotonicity of the individual size bounds is assumed.

One fixed stdout-only run,48-second deadline,20GiB free-space guard,
estimated peak below128MiB. No generated files, bytecode, or parameter scan.
"""

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="e6c058696a5383b97df74c65a1f5a18866854fe010c202da075f950fa6c2f69d"


def load_engine(deadline):
    path=Path(__file__).with_name("cheap_lane_quiet_window_cover.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen quiet-window engine hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_quiet_window_cover",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.DEADLINE=deadline
    return module


def profiles(size):
    half=size//2
    return ((half,0),(half-1,2)) if size%2==0 else ((half,1),(half,1))


def checks(engine):
    normalization=cyclic=0
    # Rank zero removes the y marker: exact coefficients must recover the
    # multinomial normalization, including every singleton orientation.
    for slots in (1,2,3):
        p=engine.quiet_half(slots,0,9,2)
        for count in (1,2,3):
            power=engine.positive_power(p,count,9,2)
            total=slots*count+4
            for f in range(10):
                for s in range(3):
                    got=engine.extract_with_free(power,f,s,4)
                    expected=(engine.DEN**count*math.comb(total,f)*math.comb(total-f,s)
                              if f+s<=total else 0)
                    assert got==expected
                    normalization+=1
    # Direct cyclic matching counts check even and odd profile formulas.
    for length in (8,10,12):
        for size in range(1,length):
            for start in range(length):
                chosen={(start+i)%length for i in range(size)}
                actual=[]
                for shift in (0,1):
                    occupancies=[int((i+shift)%length in chosen)
                                 +int((i+shift+1)%length in chosen)
                                 for i in range(0,length,2)]
                    actual.append((occupancies.count(2),occupancies.count(1)))
                assert sorted(actual)==sorted(profiles(size))
                cyclic+=1
    engine.guard()
    return {"alpha_zero_coefficient_checks":normalization,"cyclic_support_profile_checks":cyclic}


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed exact band check accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic()
    engine=load_engine(started+48)
    tests=checks(engine)
    poly=engine.positive_power(engine.quiet_half(),166)
    common_denominator=engine.DEN**166
    regions={}
    for size in range(7,257):
        for f,s in profiles(size):
            if (f,s) not in regions:
                regions[f,s]=(engine.extract_with_free(poly,f,s,1440),
                    common_denominator*math.comb(4096,f)*math.comb(4096-f,s))
        engine.guard()
    constant=2*8192*85**256
    margins={}
    for size in range(7,257):
        (a,b),(c,d)=(regions[p] for p in profiles(size))
        numerator=(a*c)**128*constant*4**(128*size)*256**(2*size)
        denominator=(b*d)**128*3**(128*size)*255**(2*size)
        margins[size]=engine.margin(numerator,denominator)
        engine.guard()
    # Positive dyadic sums give exact union bounds. Clip each size bound
    # at1 first; this affects no successful40-bit tail.
    ceiling_exponent=max(0,*margins.values())
    denominator=1<<ceiling_exponent
    total=0; tails={}; first=None
    for size in range(256,6,-1):
        total+=1<<(ceiling_exponent-max(0,margins[size]))
        tails[size]=(total,engine.margin(total,denominator))
        if total*(1<<40)<=denominator:
            first=size
    selected={7,16,32,64,96,128,160,192,224,256}
    if first is not None:
        selected.update((first-1,first,first+1))
    selected=sorted(r for r in selected if 7<=r<=256)
    table_bytes=json.dumps(sorted(margins.items()),separators=(",",":")).encode()
    engine.guard()
    result={
        "scope":"exact consecutive-support quiet-window size-band bound",
        "whole_code_certificate":False,"checked_every_size":[7,256],
        "checked_size_count":len(margins),"selected_quiet_packets_per_region":83,
        "window_packets":44,"window_choices_per_region":85,
        "source_intervals_per_size":8192,"rank_mark":"3/4",
        "source_factor_upper":"2*(4/3)^(128*r)*(256/255)^(2*r)",
        "target_margin_bits":40,"smallest_certified_dyadic_tail_start":first,
        "certified_tail":None if first is None else [first,256],
        "certified_tail_size_count":0 if first is None else 257-first,
        "certified_tail_union_margin_bits":None if first is None else tails[first][1],
        "preceding_tail_union_margin_bits":None if first in (None,7) else tails[first-1][1],
        "selected_per_size_integer_margins":{str(r):margins[r] for r in selected},
        "selected_dyadic_tail_margins":{str(r):tails[r][1] for r in selected},
        "all_size_margin_table_sha256":hashlib.sha256(table_bytes).hexdigest(),
        "tests":tests,"engine_sha256":PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds":time.monotonic()-started,"generated_bytes":0,
        "estimated_peak_bytes":128*(1<<20)
    }
    if first is not None:
        n=tails[first][0]
        result["tail_dyadic_numerator_sha256"]=hashlib.sha256(n.to_bytes((n.bit_length()+7)//8,"big")).hexdigest()
        result["tail_dyadic_denominator_exponent"]=ceiling_exponent
        result["tail_dyadic_numerator_bits"]=n.bit_length()
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
