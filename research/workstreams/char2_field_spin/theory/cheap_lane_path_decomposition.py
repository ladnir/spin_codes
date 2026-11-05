"""FLOAT decomposition of one regional16-to16 majorant entry.

Keep the frozen raw caps and tile4 route-Jensen rule at r256,z.85,alpha.523.
Compare three positive path classes: all packet-boundary masks are full16;
all masks are nonempty; and the unrestricted operator. The incoming and
outgoing regional masks are full16 in every case. Restrictions are imposed
BEFORE raw tile multiplication and the fractional power. This diagnostic
is an oracle restriction of a majorant, not an upper bound on the full code.
Because Jensen is nonlinear, differences between the final logarithms are
comparisons of restricted operators, not probabilities of complementary
events. The repeated fully paired marginal is not a feasible whole-chain
source class.

The restricted cases have only a nonquiet channel. Their destination
multiplicities are1 for full16 and C(16,j) for j=1..16 otherwise. No
incoming-state multiplicity is added. The unrestricted case retains the
frozen positive quiet/nonquiet channels. Regional input degree128 is
extracted after128 regular packets; no terminal is used in this comparison.

For orientation, report the typical full-state raw cap and the exact
output-only positive zero-subset moment after excluding four outside-input
zeros in an active RS block. For block occupancy ell this moment is
g^ell sum_{v<=min(3,16-ell)} C(16-ell,v)z^(16-ell-v)(1-z)^v q^-v,
where g=z+(1-z)/256. All requested-zero ranks here are exact for full H.
This local comparison does not add a new cap to any frozen operator.

One bounded invocation:48seconds,20GiB free,estimated peak<256MiB.
Run python -B. Stdout only; no generated files or frozen edits.
"""

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
TILES_SHA256 = "a96bcfc95c927501264f9787692e509b03f64a8f968570bc0ac631f581d8f34e"


def load_tiles():
    path = Path(__file__).with_name("cheap_lane_mask_tiles.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != TILES_SHA256:
        raise RuntimeError("frozen mask-tiles helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_mask_tiles", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def restricted_packet(m, states, degrees, weights, np, middle):
    sub = m[:, states, :][:, :, states]
    dense = np.zeros((33, len(states), len(states)))
    np.add.at(dense, degrees, weights[:, None, None]*sub)
    return np.full((33, len(states)), float("-inf")), middle.to_log(dense)


def regional(poly, destinations, args, tiles, paired):
    for _ in range(2):
        poly = paired.positive_convolution(poly, poly,128)
    poly = tiles.jensen(poly,128,destinations,args.alpha,paired.np)
    for _ in range(5):
        poly = paired.positive_convolution(poly, poly,128)
    q,n = poly
    denominator = math.log(math.comb(4096,128))
    return float((n[128,-1,-1]-denominator)/paired.LN2)


def full_state_moment(ell,z):
    n =16-ell
    return (z+(1-z)/256)**ell*sum(math.comb(n,v)*z**(n-v)*(1-z)**v*256**(-v)
                                  for v in range(min(3,n)+1))


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed diagnostic accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start=time.monotonic()
    tiles=load_tiles(); grouped=tiles.load_grouped(); singleton=grouped.load_singleton()
    paired=singleton.load_block().load_paired(start+48)
    tests={**tiles.toy_test(grouped,paired),**grouped.toy_test(paired)}
    args=SimpleNamespace(r=256,z=.85,alpha=.523,activity_mark=1.0)
    filtered,middle=paired.imports(); np=paired.np
    m,t,degrees,weights=grouped.raw_tables(args,singleton,paired,filtered)
    halves=[(a,b) for a in range(17) for b in range(a,17)]
    typical=[]
    for half in ((0,0),(0,1),(1,1),(0,2)):
        index=halves.index(half)
        moment=math.prod(full_state_moment(l,args.z)**2 for l in half)
        raw=float(m[index,16,16])
        typical.append({"half_counts":half,"raw_cap":raw,
            "filtered_output_only_exact_moment":moment,
            "raw_vs_z64_log2":math.log2(raw/args.z**64),
            "raw_vs_exact_filtered_moment_log2":math.log2(raw/moment)})
    results=[]
    for label,states in (("A_all_full16",[16]),("B_all_nonempty",list(range(1,17))),
                         ("C_unrestricted",list(range(17)))):
        paired.guard(); stage=time.monotonic()
        if label.startswith("C"):
            poly,_=tiles.raw_polynomials(m,t,degrees,weights,np,middle)
        else:
            poly=restricted_packet(m,states,degrees,weights,np,middle)
        destinations=[math.comb(16,j) for j in states]
        value=regional(poly,destinations,args,tiles,paired)
        results.append({"operator":label,"regional_nonquiet16to16_log2":value,
                        "seconds":time.monotonic()-stage})
    print(json.dumps({"scope":"FLOAT restricted-majorant one-region diagnostic",
        "whole_code_certificate":False,"actual_support_class_bound":False,
        "oracle_restrictions_are_not_code_upper_bounds":True,
        "parameters":vars(args),"tile_packets":4,
        "ideal_z64_fullstate_regional_log2":args.alpha*8192*math.log2(args.z),
        "typical_local_fullstate":typical,"results":results,
        "B_minus_A_log2":results[1]["regional_nonquiet16to16_log2"]-results[0]["regional_nonquiet16to16_log2"],
        "C_minus_B_log2":results[2]["regional_nonquiet16to16_log2"]-results[1]["regional_nonquiet16to16_log2"],
        "tests":tests,"estimated_peak_bytes":64*(1<<20),
        "seconds":time.monotonic()-start,"generated_bytes":0,
        "tiles_source_sha256":TILES_SHA256},indent=2))


if __name__=="__main__":
    main()
