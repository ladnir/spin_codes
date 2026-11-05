"""FLOAT fixed-window diagnosis of route concentration versus source counting.

Fix256 consecutive source groups, aligned with the even matching. They
occupy128 even pairs and127 odd pairs plus2 odd singletons. Require every
occupied routing pair to land in the first44 packets of each region.
There are1408 such pair slots out of4096. The exact route probability is

    [C(1408,128)/C(4096,128)]^128
    [C(1408,129)/C(4096,129)]^128.

All possible inputs supported in this group set then vanish outside the
windows. Impose zero state after each window: at most256*16=4096 F-linear
constraints on32768 source dimensions. A nonzero kernel vector exists,
and its output is supported in44/128<35% of all packets. This is a genuine
bad-code event, but an exceptionally rare one, not a failure of the target
probability guarantee. The vector need not use all256 groups. At most32
group projections can vanish identically on the kernel; since |F|>256,
a vector avoids every other zero projection and uses224..256 groups.

For any proper256-group set, occupied-pair counts in the two matchings are
at least128 each and cannot both equal128. Because P(k)=C(1408,k)/C(4096,k)
decreases with k, the displayed probability also upper-bounds every such
set's fixed-window event. Multiplying by C(8192,256) therefore covers their
union. Multiplying further by85^256 covers all length44 contiguous packet
windows independently chosen in each region. This is only the specified
route-concentration class, not all bad codes.

Separately, calculate a fixed-path contribution to the existing scalar
majorant, with its artificial repeated fully-paired marginal. The window
path is0->16, then42 steps16->16, then16->0. All84 outside packets have
empty input and zero state. Positive raw cap polynomials A,L,E give
[x^128] A L^42 E / C(4096,128), after powering each fixed-path packet cap
by alpha. Route counts stay outside alpha. This is a positive term bounded
by the unrestricted tile4 majorant; it is NOT the true event probability.
Its repeated-regional score deliberately retains the old source prefactor
to expose overcounting. Do not conflate that artificial score with the real
contiguous-support route probability above.

Run python -B. One48-second stdout-only invocation;20GiB free;estimated
peak<256MiB; no generated files, frozen edits, or parameter scan.
"""

from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode=True
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
os.environ.setdefault("OMP_NUM_THREADS","1")
GROUPED_SHA256="4526522c47c9920482f8990bb104243eada12318774ee2d0dbb68fda3704ab36"


def load_grouped():
    path=Path(__file__).with_name("cheap_lane_grouped_masks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=GROUPED_SHA256:
        raise RuntimeError("frozen grouped helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_grouped_masks",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def log_convolve(a,b,limit,np,guard):
    result=np.full(min(limit,len(a)+len(b)-2)+1,float("-inf"))
    for d in range(len(result)):
        guard()
        lo,hi=max(0,d-len(b)+1),min(d,len(a)-1)
        result[d]=np.logaddexp.reduce(a[lo:hi+1]+b[d-hi:d-lo+1][::-1])
    return result


def log_power(a,power,limit,np,guard):
    result=np.array([0.0])
    while power:
        if power&1:
            result=log_convolve(result,a,limit,np,guard)
        power//=2
        if power:
            a=log_convolve(a,a,limit,np,guard)
    return result


def exact_product(a,b,limit):
    out=[Fraction(0)]*(min(limit,len(a)+len(b)-2)+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<len(out):
                out[i+j]+=x*y
    return out


def tests(paired):
    np=paired.np
    p=[Fraction(1,2),Fraction(2,3),Fraction(3,7)]
    checked=0
    for power in (1,2,4):
        exact=[Fraction(1)]
        for _ in range(power):
            exact=exact_product(exact,p,5)
        log=log_power(np.log(np.array(p,dtype=float)),power,5,np,paired.guard)
        for x,y in zip(exact,log):
            assert math.isclose(math.log(x),float(y),rel_tol=1e-13,abs_tol=1e-13)
            checked+=1
    matching_cases=0
    for r in (2,4,6):
        for subset in itertools.combinations(range(8),r):
            active=set(subset)
            counts=[]
            for shift in (0,1):
                counts.append(sum(bool(active.intersection(((i+shift)%8,(i+shift+1)%8)))
                                  for i in range(0,8,2)))
            assert min(counts)>=r//2 and sum(counts)>=r+1
            for capacity in range(4):
                def probability(k):
                    return Fraction(math.comb(capacity,k),math.comb(4,k)) if k<=capacity else Fraction(0)
                assert probability(counts[0])*probability(counts[1])<=probability(r//2)*probability(r//2+1)
                matching_cases+=1
    return {"exact_scalar_convolution_checks":checked,
            "exact_small_cycle_window_envelopes":matching_cases}


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed diagnostic accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start=time.monotonic()
    grouped=load_grouped(); singleton=grouped.load_singleton()
    paired=singleton.load_block().load_paired(start+48)
    checks=tests(paired)
    filtered,middle=paired.imports(); np=paired.np
    args=SimpleNamespace(r=256,z=.85,alpha=.523,activity_mark=1.0)
    m,_,degrees,weights=grouped.raw_tables(args,singleton,paired,filtered)
    polynomials=[]
    for h,j in ((0,16),(16,16),(16,0)):
        coefficient=np.zeros(33)
        np.add.at(coefficient,degrees,weights*np.power(m[:,h,j],args.alpha))
        polynomials.append(middle.to_log(coefficient))
    a,l,e=polynomials
    live=log_power(l,42,128,np,paired.guard)
    path=log_convolve(log_convolve(a,live,128,np,paired.guard),e,128,np,paired.guard)
    log_region=float(path[128]-math.log(math.comb(4096,128)))/paired.LN2
    eta=256**2/(255**2*(257-args.r))
    logbeta=160*math.log2(256)-32*math.log2(255)
    limit=35*(1<<21)//100
    source=math.log2(math.comb(8192,args.r))+args.alpha*(math.log2(eta)+args.r*logbeta-limit*math.log2(args.z))
    logp128=math.log2(math.comb(1408,128))-math.log2(math.comb(4096,128))
    logp129=math.log2(math.comb(1408,129))-math.log2(math.comb(4096,129))
    fixed=128*(logp128+logp129)
    source_union=fixed+math.log2(math.comb(8192,256))
    all_windows=source_union+256*math.log2(85)
    print(json.dumps({"scope":"FLOAT fixed-window mechanism diagnostic",
        "whole_code_certificate":False,
        "real_route_event":{"profile_even":{"full":128,"single":0},
            "profile_odd":{"full":127,"single":2},"window_packets":44,
            "fixed_contiguous_support_log2_probability":fixed,
            "all_256_group_sets_fixed_windows_log2_union_upper":source_union,
            "all_256_group_sets_all_44_windows_log2_union_upper":all_windows,
            "source_dimension":32768,"boundary_constraints_at_most":4096,
            "kernel_dimension_at_least":28672,"some_witness_active_group_range":[224,256],
            "output_weight_at_most":44*64*256,"distance_threshold":limit},
        "artificial_fixed_path_majorant":{"matching_profile_repeated_full128":True,
            "regional_log2":log_region,"source_log2_prefactor":source,
            "repeated_regional_log2_score":source+256*log_region,
            "regular_window_transitions":{"activation":1,"full_live":42,"extinction":1},
            "outside_packets_per_region":84},
        "tests":checks,"seconds":time.monotonic()-start,"generated_bytes":0,
        "estimated_peak_bytes":32*(1<<20),"grouped_source_sha256":GROUPED_SHA256},indent=2))


if __name__=="__main__":
    main()
