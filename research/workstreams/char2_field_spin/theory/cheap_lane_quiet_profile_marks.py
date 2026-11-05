"""Three fixed rank marks for the actual size256 quiet-window profile union.

This FLOAT diagnostic imports the pinned quiet-profile screen without
modifying it. The route, window class, positive polynomial, mode smoothing,
and actual joint source-profile counts are unchanged. The three marks are
y=3/4,2/3,3/5; alpha=-log2(y)/8. A complete profile bound uses one common
mark for both matchings. The minimum is taken only between these complete
profile bounds, before summing all16640 ordered source profiles.

The sharper quiet-only source interface applies to the E-defined quiet
constraint kernel. At most256 source-group projections are required to be
nonzero. A field-extension witness with exact group support therefore
implies an E-scalar witness with exact support. The number of such scalar
witnesses is either zero or at least255, by scalar multiplication. Thus
the existence indicator is at most n_T/255. Applying the summed outer
density and the fractional bound gives the exact source log2 factor
  -32768*log2(y)+alpha*(8192*log2(256/255)-log2(255)).
No eta comparison factor is required for this quiet-only interface.

Outputs are FLOAT diagnostics, never a certificate. A negative union
requires independent directed/exact reevaluation. One fixed48-second run,
20GiB disk guard, estimated peak below64MiB, stdout only; no generated data.
"""

import hashlib
import heapq
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="28083dbe863b0408c2f5b7c22cbc09302f444bdb4682d1d2bb9ce1af7072e151"
MARKS=(("3/4",3/4),("2/3",2/3),("3/5",3/5))


def load_screen(deadline):
    path=Path(__file__).with_name("cheap_lane_quiet_profile_screen.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen profile screen hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_quiet_profile_screen",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.DEADLINE=deadline
    return module


def keep_worst(heap,item,limit=6):
    if len(heap)<limit:
        heapq.heappush(heap,item)
    elif item>heap[0]:
        heapq.heapreplace(heap,item)


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed three-mark screen accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); screen=load_screen(started+48)
    counts=screen.load_counts()
    checks=screen.checks(screen.terms(.75))
    data=[]
    for name,y in MARKS:
        rows=screen.terms(y); regional=[]
        for full in range(129):
            row=screen.optimize(rows,full)
            row["multinomial_mode_log2"]=screen.multinomial_mode_log(row["u"],row["v"])/screen.LN2
            row["regional_smoothed_log2"]=row["regional_plain_log2"]+row["multinomial_mode_log2"]
            regional.append(row); screen.guard()
        alpha=-math.log2(y)/8
        source=-32768*math.log2(y)+alpha*(8192*math.log2(256/255)-math.log2(255))
        data.append({"mark":name,"alpha":alpha,"source_log2":source,
                     "common":source+256*math.log2(85),"regional":regional})
    per_mark={name:{"log2_sum":-math.inf,"positive_profiles":0,"negative_profiles":0,
                    "selected_by_complete_minimum":0,"worst":[]} for name,_ in MARKS}
    total=-math.inf; worst=[]; profile_count=0; selected_positive=selected_negative=0
    selected_runs={c:{"count":0,"minimum":math.inf,"maximum":-math.inf}
                   for c in (1,4,16,64,128,196,256)}
    for f0 in range(129):
        for f1 in range(129):
            if f0+f1>255:
                continue
            count=counts.profile_count_upper(4096,256,f0,f1)
            if not count:
                continue
            c=256-f0-f1; profile_count+=1; candidates=[]
            for mark in data:
                score=mark["common"]+math.log2(count)+128*(mark["regional"][f0]["regional_smoothed_log2"]
                                                         +mark["regional"][f1]["regional_smoothed_log2"])
                name=mark["mark"]; result=per_mark[name]
                result["log2_sum"]=screen.logadd(result["log2_sum"],score)
                result["positive_profiles" if score>=0 else "negative_profiles"]+=1
                keep_worst(result["worst"],(score,f0,f1,c))
                candidates.append((score,name))
            score,name=min(candidates)
            per_mark[name]["selected_by_complete_minimum"]+=1
            total=screen.logadd(total,score)
            selected_positive+=score>=0; selected_negative+=score<0
            keep_worst(worst,(score,f0,f1,c,name))
            if c in selected_runs:
                row=selected_runs[c]; row["count"]+=1
                row["minimum"]=min(row["minimum"],score); row["maximum"]=max(row["maximum"],score)
        screen.guard()
    for value in per_mark.values():
        value["worst"]=[{"log2_bound":score,"full0":f0,"full1":f1,"runs":c}
                        for score,f0,f1,c in sorted(value["worst"],reverse=True)]
    witness_indices={name:set() for name,_ in MARKS}
    for _,f0,f1,_,name in worst:
        witness_indices[name].update((f0,f1))
    print(json.dumps({"scope":"FLOAT size256 actual-profile quiet-window three-mark union",
        "whole_code_certificate":False,"directed":False,"ordered_profile_count":profile_count,
        "log2_sum_of_complete_profile_minima":total,
        "selected_positive_profiles":selected_positive,"selected_negative_profiles":selected_negative,
        "worst_selected_profiles":[{"log2_bound":score,"full0":f0,"full1":f1,"runs":c,"mark":name}
                                   for score,f0,f1,c,name in sorted(worst,reverse=True)],
        "per_mark":per_mark,"selected_run_counts":selected_runs,
        "mark_parameters":[{"mark":m["mark"],"alpha":m["alpha"],"source_log2":m["source_log2"],
            "optimizer_converged_profiles":sum(r["converged"] for r in m["regional"]),
            "maximum_gradient":max(r["maximum_gradient"] for r in m["regional"]),
            "worst_profile_regional_witnesses":[m["regional"][i] for i in sorted(witness_indices[m["mark"]])]}
            for m in data],
        "tests":checks,"source_interface":"quiet exact-support scalar projectivization; no eta",
        "screen_sha256":PIN,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "seconds":time.monotonic()-started,"generated_bytes":0,"estimated_peak_bytes":64*(1<<20)
    },indent=2))


if __name__=="__main__":
    main()
