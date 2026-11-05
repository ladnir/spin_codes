"""Controlled raw-pooling comparison at16,64,128 packets; FLOAT only.

Import the frozen degree-conditioned tile helper and retain its raw caps,
actual consecutive256-group source profiles, z=.85, alpha=.05, quiet mark,
floor mark1, and OLD route witnesses. No parameter is optimized here.
Raw positive log polynomials reach p packets by doubling P and P^(p-1)T.
At p, normalize each coefficient by C(32p,f)C(32p-f,s)*D_j before alpha;
restore route and endpoint multiplicities afterward. The output-only
terminal has the same32p route slots but only p-1 updates and no final D.

Evaluate every checkpoint at the SAME OLD route witnesses. Regional
assembly uses128/p tiles, with the last tile containing the terminal.
For p128, a region is a single tile and its terminal is P127T. Final
assembly retains both actual alternating orientations and32767 updates.

Report three channels. The full17-state channel is the majorant. The
all-packet0/16 channel restricts every raw incoming/outgoing state before
any polynomial power. The boundary0/16 channel restricts only tile
endpoints after raw pooling, leaving intermediate states unrestricted.
The two restrictions are positive submajorants, NOT upper bounds for the
full code or the full source class. Pool-size improvements combine fewer
boundary mask counts, larger Jensen groups, and stronger degree cutoffs.
At p16, evaluating eight tiles still admits impossible regional degrees;
at p128 those total degrees have already been truncated. Consequently the
comparison is not a pure isolation of boundary-mask entropy.

One90-second guarded job,20GiB free,128MiB estimated peak, stdout only.
No generated files, caches, frozen edits, scans, optimizers, or retries.
Progress is flushed at every raw doubling so timeout evidence survives.
"""

import ctypes
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
PIN="264a216fd5beb0c2b2d18826e5fc3f7cd4638561ff2aa37d9862382e36bbd0d2"
LN2=math.log(2)
DEADLINE=0.
np=None
lse=None


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("90-second deadline exhausted; no retry")


def load():
    path=Path(__file__).with_name("cheap_lane_degree_conditioned_tile.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen degree-conditioned helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_pooling_degree",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.np=np; module.lse=lse; module.DEADLINE=DEADLINE
    return module,module.load()


def peak_working_set():
    if os.name!="nt":
        return None
    class Counters(ctypes.Structure):
        _fields_=[("cb",ctypes.c_ulong),("faults",ctypes.c_ulong)]+[
            (name,ctypes.c_size_t) for name in ("peak","working","paged_peak","paged",
                "nonpaged_peak","nonpaged","pagefile","pagefile_peak")]
    value=Counters(); value.cb=ctypes.sizeof(value)
    kernel=ctypes.WinDLL("kernel32",use_last_error=True)
    kernel.GetCurrentProcess.restype=ctypes.c_void_p
    psapi=ctypes.WinDLL("psapi",use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.POINTER(Counters),ctypes.c_ulong]
    if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(value),value.cb):
        return None
    return int(value.peak)


def tile_value(poly,slots,logu,logw,logU):
    pieces=[]; full=np.arange(poly.shape[1],dtype=float)
    for single in range(len(poly)) if math.isfinite(logw) else range(1):
        pieces.append(poly[single]+(full*logu+(single*logw if single else 0))[:,None,None])
    return lse(np.concatenate(pieces,axis=0),axis=0)-slots*logU


def score(raw,terminal,packets,channel,degree,pilot):
    guard(); slots=32*packets
    d=[1,1] if raw.shape[2]==2 else [math.comb(16,j) for j in range(17)]
    a=degree.normalize(raw,slots,degree.ALPHA,d)
    b=degree.normalize(terminal,slots,degree.ALPHA,[1])
    if channel=="boundary0_16":
        a=a[:,:,[0,16],:][:,:,:,[0,16]]
        b=b[:,:,[0,16],:]
    regions=[]; tails=[]; witnesses=[]
    logs=np.log(np.asarray(degree.OLD))
    for i,(f,s) in enumerate(((128,0),(127,2))):
        logu=float(logs[i]); logw=-math.inf if i==0 else float(logs[2])
        u=math.exp(logu); w=0. if i==0 else math.exp(logw); logU=math.log1p(u+w)
        tile=tile_value(a,slots,logu,logw,logU)
        tail=tile_value(b,slots,logu,logw,logU)
        reg,end=pilot.double_pair(tile,tail,7-(packets.bit_length()-1))
        regions.append(reg); tails.append(end)
        h=4096*logU-f*logu-(s*logw if s else 0)-math.log(math.comb(4096,f)*math.comb(4096-f,s))
        last=reg.shape[0]-1
        witnesses.append({"full":f,"single":s,"H_bits":h/LN2,
            "regular00_bits":float(reg[0,0]/LN2),"regular16to16_bits":float(reg[last,last]/LN2),
            "terminal0_bits":float(end[0,0]/LN2)})
    paths,orientations=pilot.chain(regions,tails)
    source=degree.ALPHA*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
    output=-degree.ALPHA*pilot.LIMIT*math.log(degree.Z)
    correction=-degree.ALPHA*24575*math.log(degree.QUIET)
    routing=128*sum(row["H_bits"] for row in witnesses)*LN2
    value=(source+output+correction+routing+math.log(4096)+paths)/LN2
    return {"pool_packets":packets,"channel":channel,"full_class_upper_bound":channel=="all17",
        "log2_value":value,"components_bits":{"source":source/LN2,"output":output/LN2,
            "quiet_correction":correction/LN2,"routing_normalizer":routing/LN2,
            "sets_per_orientation":12.,"finite_two_orientation_sum":paths/LN2},
        "orientation_chain_bits":[x/LN2 for x in orientations],"regional":witnesses}


def main():
    global DEADLINE,np,lse
    if len(sys.argv)!=1:
        raise ValueError("fixed pooling comparison accepts no scans")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); DEADLINE=started+90
    import numpy as numpy
    from scipy.special import logsumexp
    np=numpy; lse=logsumexp
    degree,pilot=load(); checked=degree.tests(pilot)
    singleton,filtered=pilot.imports()
    shapes,index,degrees,coeff,records,count_checks=pilot.geometry()
    m,t,refined=pilot.raw_tables(shapes,index,records,degree.Z,singleton,filtered)
    a,b=degree.raw_polynomials(m,t,degrees,coeff,pilot)
    small_a=a[:,:,[0,16],:][:,:,:,[0,16]].copy()
    small_b=b[:,:,[0,16],:].copy()
    del m,t
    results=[]; stages=[]
    for stage in range(7):
        guard(); current=time.monotonic()
        b,a=degree.convolution(a,b),degree.convolution(a,a)
        small_b,small_a=degree.convolution(small_a,small_b),degree.convolution(small_a,small_a)
        packets=2**(stage+1)
        stages.append({"packets":packets,"seconds":time.monotonic()-current,
            "total_seconds":time.monotonic()-started,"peak_working_set_bytes":peak_working_set()})
        print(json.dumps({"progress":"raw_pool_complete",**stages[-1]}),flush=True)
        if packets in (16,64,128):
            full=score(a,b,packets,"all17",degree,pilot)
            boundary=score(a,b,packets,"boundary0_16",degree,pilot)
            small=score(small_a,small_b,packets,"all_packet0_16",degree,pilot)
            assert small["log2_value"]<=boundary["log2_value"]+1e-7
            assert boundary["log2_value"]<=full["log2_value"]+1e-7
            if packets==16:
                assert abs(full["log2_value"]-36376.29710010885)<1e-6
            results.extend((full,boundary,small))
            print(json.dumps({"progress":"checkpoint","pool_packets":packets,
                "log2_values":{row["channel"]:row["log2_value"] for row in (full,boundary,small)},
                "seconds":time.monotonic()-started}),flush=True)
    guard()
    print(json.dumps({"scope":"FLOAT fixed-witness pooling comparison on actual consecutive256 supports; restricted channels are submajorants",
        "whole_code_certificate":False,"directed":False,
        "interpretation":"pool-size changes combine boundary counts, Jensen grouping, and stronger degree truncation",
        "parameters":{"z":degree.Z,"alpha":degree.ALPHA,"quiet_mark":degree.QUIET,"floor_mark":1,
                      "route_witness":degree.OLD,"quiet_count_upper":24575},
        "results":results,"raw_stages":stages,"tests":{**checked,"exact_packet_route_normalizers":count_checks,
            "positive_channel_inclusions":6,"frozen_pool16_reproduction":"passed"},
        "terminal_updates_by_pool":{"16":15,"64":63,"128":127},"charged_regular_steps":32767,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"degree_helper_pin":PIN,
        "estimated_peak_bytes":128*(1<<20),"peak_working_set_bytes":peak_working_set(),
        "seconds":time.monotonic()-started,"generated_bytes":0},indent=2))


if __name__=="__main__":
    main()
